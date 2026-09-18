#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Federated GWAS: server side.

Sites run their own GWAS and standardise the result to GWAMA format with
gwas_cohort_qc_with_gwama.R, so the server never converts anything and never
needs to know which tool produced the numbers. It validates what arrives,
writes it to disk, and hands the file list to run_gwama.sh.
"""

import argparse
import os
import subprocess

from nvflare.apis.fl_constant import FLContextKey
from nvflare.apis.fl_context import FLContext
from nvflare.app_common.aggregators.model_aggregator import ModelAggregator
from nvflare.app_opt.pt.recipes.fedavg import FedAvgRecipe
from nvflare.client import FLModel
from nvflare.recipe import ProdEnv, SimEnv, add_experiment_tracking

from model import DummyModel

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TOOLS_ROOT = "/home/ubuntu/tools"
WRAPPER_NAME = "run_gwama.sh"

# Every method uses the same site driver, which dispatches on PROGRAM.
METHODS = ("regenie", "saige")
SITE_SCRIPT = "local_script_start_gwas.sh"
QC_SCRIPT = "gwas_cohort_qc_with_gwama.R"

# Columns the R script writes, and the mode run_gwama.sh needs.
TRAIT_TYPES = {
    "binary": ("or", ["MARKERNAME", "EA", "NEA", "EAF", "N",
                      "OR", "OR_95L", "OR_95U"]),
    "quantitative": ("qt", ["MARKERNAME", "EA", "NEA", "EAF", "N",
                            "BETA", "SE"]),
}


def _run_dir(fl_ctx: FLContext) -> str:
    workspace = fl_ctx.get_prop(FLContextKey.WORKSPACE_OBJECT)
    return workspace.get_run_dir(fl_ctx.get_job_id())


def _find_wrapper():
    """
    Locate run_gwama.sh at run time.

    In production the aggregator executes inside the job's custom directory on
    the FL server, not in the repo it was submitted from, so a path baked in
    at import time points at the wrong place.
    """
    candidates = [
        os.environ.get("FEDGX_GWAMA_WRAPPER"),
        os.path.join(SCRIPT_DIR, WRAPPER_NAME),
        os.path.join(os.getcwd(), WRAPPER_NAME),
        os.path.join(DEFAULT_TOOLS_ROOT, WRAPPER_NAME),
    ]
    for path in candidates:
        if path and os.path.isfile(path):
            return path
    raise FileNotFoundError(
        f"{WRAPPER_NAME} not found in {[c for c in candidates if c]}; "
        "set FEDGX_GWAMA_WRAPPER"
    )


class GWASMetaAggregator(ModelAggregator):
    """Collect GWAMA-format summary statistics and meta-analyse them."""

    def __init__(self, trait_type="binary", model="fixed",
                 tools_root=DEFAULT_TOOLS_ROOT, output_folder="server_results"):
        super().__init__()
        self.trait_type = trait_type
        self.model = model
        self.tools_root = tools_root
        self.output_folder = output_folder

        self.received_params_type = None
        self.passthrough_params = None
        self.output_dir = None
        self.accepted = []          # (dataset_id, path)
        self.rejected = []          # (site_name, reason)

    def _check(self, path):
        """Return None if the file is usable, else why it isn't."""
        _, required = TRAIT_TYPES[self.trait_type]
        try:
            with open(path) as f:
                header = f.readline().rstrip("\n").split("\t")
                n_rows = sum(1 for _ in f)
        except OSError as e:
            return f"unreadable: {e}"

        missing = [c for c in required if c not in header]
        if missing:
            return f"missing columns {missing} for trait type {self.trait_type}"
        if n_rows == 0:
            return "no data rows"
        return None

    def accept_model(self, model: FLModel):
        if self.received_params_type is None:
            self.received_params_type = model.params_type
            # Kept only so the aggregated model has a valid state dict for the
            # persistor; nothing here is trained.
            self.passthrough_params = model.params

        if self.output_dir is None:
            self.output_dir = os.path.join(_run_dir(self.fl_ctx),
                                           self.output_folder)
            os.makedirs(self.output_dir, exist_ok=True)
            print(f"Output directory: {self.output_dir}")

        meta = model.meta or {}
        site = meta.get("site_name", "unknown_site")
        dataset_id = str(meta.get("dataset_id", "unknown_id"))

        def reject(reason):
            print(f"REJECTED {site}: {reason}")
            self.rejected.append((site, reason))

        if meta.get("success") is False:
            return reject(meta.get("error_message", "client reported failure"))

        content = meta.get("results_file", "")
        if not content:
            return reject("no results_file content")

        site_trait = meta.get("trait_type")
        if site_trait and site_trait != self.trait_type:
            return reject(f"ran trait_type '{site_trait}', "
                          f"server expects '{self.trait_type}'")

        path = os.path.join(self.output_dir, f"site{dataset_id}_{site}_gwama.txt")
        with open(path, "w") as f:
            f.write(content)

        qc = meta.get("qc_summary", "")
        if qc:
            with open(os.path.join(self.output_dir,
                                   f"site{dataset_id}_{site}_qc.txt"), "w") as f:
                f.write(qc)

        problem = self._check(path)
        if problem:
            return reject(problem)

        self.accepted.append((dataset_id, path))
        print(f"Accepted {site} (site{dataset_id}): {path}")

    def aggregate_model(self) -> FLModel:
        for site, reason in self.rejected:
            print(f"Excluded {site}: {reason}")

        if len(self.accepted) < 2:
            raise RuntimeError(
                f"only {len(self.accepted)} usable site file(s); "
                "a meta-analysis needs at least 2"
            )
        wrapper = _find_wrapper()
        mode, _ = TRAIT_TYPES[self.trait_type]
        output_root = os.path.join(self.output_dir, "gwama")

        # Sorted: GWAMA takes the reference allele from the first file in the
        # list, so the order must not depend on which site replied first.
        files = [path for _, path in sorted(self.accepted)]
        cmd = ["bash", wrapper, mode, "--model", self.model,
               output_root, *files]

        env = os.environ.copy()
        env.setdefault("FEDGX_TOOLS_ROOT", self.tools_root)

        print(f"\nMeta-analysis over {len(files)} sites")
        print("Command: " + " ".join(cmd))

        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
        if result.returncode != 0:
            raise RuntimeError(f"run_gwama.sh exited {result.returncode}")

        return FLModel(
            params=self.passthrough_params,
            params_type=self.received_params_type,
            meta={
                "META_ANALYSIS_COMPLETED": True,
                "N_SITES": len(files),
                "N_REJECTED": len(self.rejected),
            },
        )

    def reset_stats(self):
        self.accepted = []
        self.rejected = []
        self.received_params_type = None
        self.passthrough_params = None


def define_parser():
    p = argparse.ArgumentParser()
    p.add_argument("--n_clients", type=int, default=3)
    p.add_argument("--num_rounds", type=int, default=1)
    p.add_argument("--env", default="sim", choices=["sim", "prod"],
                   help="default: %(default)s")
    p.add_argument("--method", default="regenie", choices=METHODS,
                   help="GWAS tool the clients run (default: %(default)s)")
    p.add_argument("--trait_type", default="binary", choices=sorted(TRAIT_TYPES),
                   help="passed to the site R script (default: %(default)s)")
    p.add_argument("--model", default="fixed",
                   choices=["fixed", "random", "both"],
                   help="meta-analysis model (default: %(default)s)")
    p.add_argument("--tools_root",
                   default=os.environ.get("FEDGX_TOOLS_ROOT", DEFAULT_TOOLS_ROOT),
                   help="where GWAMA lives on the server (default: %(default)s)")
    p.add_argument("--startup_kit", default=os.environ.get("FEDGX_STARTUP_KIT"),
                   help="prod only (or set FEDGX_STARTUP_KIT)")
    p.add_argument("--username", default=os.environ.get("FEDGX_USERNAME"),
                   help="prod only (or set FEDGX_USERNAME)")

    args = p.parse_args()
    if args.env == "prod" and not (args.startup_kit and args.username):
        p.error("--env prod requires --startup_kit and --username "
                "(or FEDGX_STARTUP_KIT / FEDGX_USERNAME)")
    return args


def main():
    args = define_parser()

    recipe = FedAvgRecipe(
        name="fed_gwas",
        model=DummyModel(),
        min_clients=args.n_clients,
        num_rounds=args.num_rounds,
        train_script="client.py",
        train_args=f"--method {args.method} --trait_type {args.trait_type}",
        aggregator=GWASMetaAggregator(trait_type=args.trait_type,
                                      model=args.model,
                                      tools_root=args.tools_root),
    )
    add_experiment_tracking(recipe, tracking_type="tensorboard")

    # Ship the site driver and the shared QC script, so every site
    # standardises with identical code.
    recipe.job.to_clients(SITE_SCRIPT)
    recipe.job.to_clients(QC_SCRIPT)

    # The aggregator runs on the FL server inside the job directory, so the
    # meta-analysis wrapper has to travel with it.
    recipe.job.to_server(WRAPPER_NAME)

    if args.env == "sim":
        print(f"Simulation environment, {args.n_clients} clients")
        env = SimEnv(num_clients=args.n_clients)
    else:
        print(f"Production environment, startup kit {args.startup_kit}")
        env = ProdEnv(startup_kit_location=args.startup_kit,
                      username=args.username, login_timeout=300)

    run = recipe.execute(env)
    print("\nJob status:", run.get_status())
    print("Results in:", run.get_result())


if __name__ == "__main__":
    main()