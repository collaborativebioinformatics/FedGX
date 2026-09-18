#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NVFLARE client for federated GWAS.

Per round the client:
  1. resolves the dataset id it was given against a site-local dataset registry
  2. resolves the GWAS tool name against a site-local tool registry
  3. validates both (files exist, binary is executable)
  4. runs local_script_start_gwas.sh, which runs the GWAS tool and then
     standardises the result to GWAMA format with gwas_cohort_qc_with_gwama.R
  5. returns the GWAMA-format text to the server

All standardisation happens here, at the site. The server receives GWAMA
format and never needs to know which tool produced it.

Both registries are site-local: nothing about paths or binaries comes from the
server, which only ever sends a dataset id and a method name.

datasets.json  (FEDGX_DATASETS, else <data_root>/datasets.json)
----------------------------------------------------------------
This file is site-local and MUST describe exactly one dataset: its own. A
registry holding several sites' entries would put every site's data layout on
every machine, so more than one entry is rejected rather than ignored.

{
  "site1": {
    "populationid": "site1",
    "bfile":  "/home/ubuntu/data/site1/site1_geno",
    "pheno":  "/home/ubuntu/data/site1/site1_pheno.pheno",
    "covar":  "/home/ubuntu/data/site1/site1_geno.covar",
    "outdir": "/home/ubuntu/data/site1/gwas_out",
    "build":  "GRCh38"
  }
}

"bfile" is a PLINK prefix without the .bed/.bim/.fam extension. "outdir" and
"build" are optional; everything else is required.

The dataset id the server sends must match the single local id, so a server
asking for someone else's dataset fails here rather than silently running on
whatever this site happens to hold.

tools.json  (FEDGX_TOOLS, else /home/ubuntu/tools/tools.json, else next to
this script)
------------------------------------------------------------------------------
Site-local, like datasets.json, but it holds only binary locations, so sites
can differ freely and listing a tool a site does not have is harmless:
resolution happens only for the method actually requested.

{
  "regenie": "/home/ubuntu/tools/regenie",
  "plink":   "/home/ubuntu/tools/plink",
  "plink2":  "/home/ubuntu/tools/plink2",
  "gcta":    "/home/ubuntu/tools/gcta64",
  "saige":   "/home/ubuntu/tools/step2_SPAtests.R",
  "GWAMA":   "/home/ubuntu/tools/GWAMA",
  "Rscript": "/usr/bin/Rscript"
}
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import traceback

import nvflare.client as flare
from nvflare.app_common.abstract.fl_model import FLModel

DEFAULT_DATA_ROOT = "/home/ubuntu/data"
DEFAULT_TOOLS_ROOT = "/home/ubuntu/tools"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GWAS_SCRIPT = os.path.join(SCRIPT_DIR, "local_script_start_gwas.sh")

REQUIRED_DATASET_KEYS = ("populationid", "bfile", "pheno", "covar")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", default="regenie",
                        help="GWAS tool to run (passed through as PROGRAM)")
    parser.add_argument("--trait_type", default="binary",
                        choices=["binary", "quantitative"])
    parser.add_argument("--data_root",
                        default=os.environ.get("FEDGX_DATA_ROOT",
                                               DEFAULT_DATA_ROOT))
    parser.add_argument("--registry", default=os.environ.get("FEDGX_DATASETS"),
                        help="dataset registry JSON "
                             "(default: <data_root>/datasets.json)")
    parser.add_argument("--tools_root",
                        default=os.environ.get("FEDGX_TOOLS_ROOT",
                                               DEFAULT_TOOLS_ROOT),
                        help="directory holding the GWAS binaries "
                             "(default: %(default)s)")
    parser.add_argument("--tools", default=os.environ.get("FEDGX_TOOLS"),
                        help="tool registry JSON "
                             "(default: <tools_root>/tools.json)")
    parser.add_argument("--timeout", type=int,
                        default=int(os.environ.get("FEDGX_GWAS_TIMEOUT",
                                                   6 * 60 * 60)),
                        help="seconds to allow the GWAS to run (default: 6h)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Registries
# ---------------------------------------------------------------------------

def _load_json(path, what):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{what} not found: {path}")
    with open(path) as f:
        try:
            return json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"{what} at {path} is not valid JSON: {e}") from e


def load_dataset_registry(args):
    """
    Load the site-local dataset registry, which must describe exactly one
    dataset: this site's own.
    """
    path = args.registry or os.path.join(args.data_root, "datasets.json")
    registry = _load_json(path, "dataset registry")

    if not isinstance(registry, dict) or not registry:
        raise ValueError(
            f"{path} must be a JSON object holding one dataset entry"
        )

    if len(registry) != 1:
        raise ValueError(
            f"{path} holds {len(registry)} dataset entries; it must hold "
            "exactly one. Each site keeps only its own dataset definition, so "
            "no site holds another site's data layout."
        )

    return registry, path


def load_tool_registry(args):
    """
    Look in the explicit location, then next to this script, then under the
    data root. A missing tool registry is not fatal: tools are then expected
    to be on PATH.
    """
    if args.tools:
        if not os.path.isfile(args.tools):
            raise FileNotFoundError(f"tool registry not found: {args.tools}")
        return _load_json(args.tools, "tool registry"), args.tools

    candidates = [
        os.path.join(args.tools_root, "tools.json"),
        os.path.join(SCRIPT_DIR, "tools.json"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return _load_json(path, "tool registry"), path

    print(f"[client] no tools.json in {candidates}; falling back to "
          f"{args.tools_root} and PATH", flush=True)
    return {}, None


def resolve_tool(tools, tools_path, name, tools_root):
    """
    Map a tool name to an executable path, in order of preference:
      1. an explicit entry in tools.json
      2. <tools_root>/<name>
      3. PATH
    """
    configured = tools.get(name)
    if configured:
        configured = os.path.expanduser(configured)
        if not os.path.isfile(configured):
            raise FileNotFoundError(
                f"'{name}' is listed in {tools_path} as {configured}, "
                "which does not exist"
            )
        if not os.access(configured, os.X_OK):
            raise PermissionError(f"{configured} is not executable")
        return configured

    in_root = os.path.join(tools_root, name)
    if os.path.isfile(in_root) and os.access(in_root, os.X_OK):
        return in_root

    found = shutil.which(name)
    if found:
        return found

    raise FileNotFoundError(
        f"'{name}' is not in {tools_path or 'tools.json'}, not at {in_root}, "
        "and not on PATH"
    )


def resolve_dataset(registry, registry_path, dataset_id, data_root):
    """
    Look up dataset_id and check the entry before anything is executed.
    """
    local_id = next(iter(registry))

    if dataset_id != local_id:
        raise KeyError(
            f"the server asked for dataset '{dataset_id}', but this site "
            f"holds '{local_id}'. Refusing to substitute."
        )

    entry = registry[local_id]

    missing_keys = [k for k in REQUIRED_DATASET_KEYS if not entry.get(k)]
    if missing_keys:
        raise ValueError(
            f"dataset '{dataset_id}' in {registry_path} is missing "
            f"required key(s): {missing_keys}"
        )

    bfile = os.path.abspath(os.path.expanduser(entry["bfile"]))
    pheno = os.path.abspath(os.path.expanduser(entry["pheno"]))
    covar = os.path.abspath(os.path.expanduser(entry["covar"]))
    outdir = os.path.abspath(os.path.expanduser(
        entry.get("outdir", os.path.join(os.path.dirname(bfile), "gwas_out"))
    ))

    # Keep resolved paths inside the configured data root unless explicitly
    # allowed otherwise. Set FEDGX_ALLOW_ANY_PATH=1 for sites that store data
    # outside it.
    if os.environ.get("FEDGX_ALLOW_ANY_PATH", "0") != "1":
        root = os.path.abspath(data_root)
        for label, path in (("bfile", bfile), ("pheno", pheno),
                            ("covar", covar)):
            if os.path.commonpath([path, root]) != root:
                raise ValueError(
                    f"dataset '{dataset_id}' {label} path {path} is outside "
                    f"the data root {root} "
                    "(set FEDGX_ALLOW_ANY_PATH=1 to permit this)"
                )

    required_files = [f"{bfile}.bed", f"{bfile}.bim", f"{bfile}.fam",
                      pheno, covar]
    missing = [p for p in required_files if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError(f"dataset '{dataset_id}' is missing: {missing}")

    return {
        "populationid": entry["populationid"],
        "bfile": bfile,
        "pheno": pheno,
        "covar": covar,
        "outdir": outdir,
        "build": entry.get("build", "GRCh38"),
    }


# ---------------------------------------------------------------------------
# Running the GWAS
# ---------------------------------------------------------------------------

def run_gwas(dataset, method, trait_type, tool_bin, rscript_bin, timeout):
    """
    Run the site GWAS script. Returns the path to the GWAMA-format output.
    """
    if not os.path.isfile(GWAS_SCRIPT):
        raise FileNotFoundError(f"GWAS script not found: {GWAS_SCRIPT}")

    env = os.environ.copy()
    env.update({
        "PROGRAM": method,
        "POPULATIONID": dataset["populationid"],
        "TRAITTYPE": trait_type,
        "BUILD": dataset["build"],
        "BFILE": dataset["bfile"],
        "PHENO_FILE": dataset["pheno"],
        "COVAR_FILE": dataset["covar"],
        "OUTDIR": dataset["outdir"],
        "TOOL_BIN": tool_bin,
        "RSCRIPT_BIN": rscript_bin,
    })

    print(f"[client] running {os.path.basename(GWAS_SCRIPT)} for "
          f"{dataset['populationid']} (method={method}, trait={trait_type})",
          flush=True)
    print(f"[client]   binary : {tool_bin}", flush=True)
    print(f"[client]   bfile  : {dataset['bfile']}", flush=True)
    print(f"[client]   pheno  : {dataset['pheno']}", flush=True)
    print(f"[client]   covar  : {dataset['covar']}", flush=True)

    proc = subprocess.run(
        ["bash", GWAS_SCRIPT],
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    # Always show the tool output; a zero exit code is not on its own evidence
    # that the GWAS produced anything usable.
    print(proc.stdout, flush=True)
    if proc.stderr:
        print(proc.stderr, file=sys.stderr, flush=True)

    if proc.returncode != 0:
        raise RuntimeError(
            f"{os.path.basename(GWAS_SCRIPT)} exited {proc.returncode}"
        )

    gwama_file = os.path.join(
        dataset["outdir"], f"{dataset['populationid']}.GWAMA.txt"
    )
    if not os.path.isfile(gwama_file) or os.path.getsize(gwama_file) == 0:
        raise RuntimeError(f"no GWAMA output produced at {gwama_file}")

    return gwama_file


def read_outputs(dataset, gwama_file):
    with open(gwama_file) as f:
        results = f.read()

    n_variants = results.count("\n") - 1
    print(f"[client] {gwama_file}: {n_variants} variants, "
          f"{len(results)} bytes", flush=True)

    qc_summary = ""
    qc_path = os.path.join(
        dataset["outdir"], f"{dataset['populationid']}.qc_summary.txt"
    )
    if os.path.isfile(qc_path):
        with open(qc_path) as f:
            qc_summary = f.read()

    return results, qc_summary, n_variants


# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    flare.init()
    site_name = flare.get_site_name()

    while flare.is_running():
        input_model = flare.receive()
        if input_model is None:
            break

        meta = input_model.meta or {}
        # A server may request an explicit dataset. Otherwise use the one and
        # only entry in this site's local registry; FLARE site names (for
        # example site-1) need not equal cohort identifiers (for example
        # site1).
        requested_dataset_id = meta.get("dataset_id")
        dataset_id = str(requested_dataset_id or site_name)
        method = meta.get("method", args.method)
        trait_type = meta.get("trait_type", args.trait_type)

        try:
            registry, registry_path = load_dataset_registry(args)
            dataset_id = str(requested_dataset_id or next(iter(registry)))
            print(f"\n[client] === {site_name} | dataset={dataset_id} | "
                  f"round={input_model.current_round} ===", flush=True)
            dataset = resolve_dataset(registry, registry_path,
                                      dataset_id, args.data_root)

            tools, tools_path = load_tool_registry(args)
            tool_bin = resolve_tool(tools, tools_path, method,
                                    args.tools_root)
            rscript_bin = resolve_tool(tools, tools_path, "Rscript",
                                       args.tools_root)

            gwama_file = run_gwas(dataset, method, trait_type,
                                  tool_bin, rscript_bin, args.timeout)
            results, qc_summary, n_variants = read_outputs(dataset, gwama_file)

            output_model = FLModel(
                # FedAvgRecipe requires a valid model state to make the
                # round-trip even though FedGX does not train it. Preserve the
                # dummy state and carry workflow status in metadata.
                params=input_model.params,
                params_type=input_model.params_type,
                meta={
                    "success": True,
                    "site_name": site_name,
                    "dataset_id": dataset_id,
                    "method": method,
                    "trait_type": trait_type,
                    "build": dataset["build"],
                    "n_variants": n_variants,
                    "results_file": results,
                    "qc_summary": qc_summary,
                },
            )

        except Exception as e:                     # noqa: BLE001
            traceback.print_exc()
            output_model = FLModel(
                params=input_model.params,
                params_type=input_model.params_type,
                meta={
                    "success": False,
                    "site_name": site_name,
                    "dataset_id": dataset_id,
                    "error_message": f"{type(e).__name__}: {e}",
                },
            )

        flare.send(output_model)


if __name__ == "__main__":
    main()
