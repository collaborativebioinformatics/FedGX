#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Importable NVFLARE aggregator and central FedGX post-processing steps."""

import os
import re
import subprocess
import sys

from nvflare.apis.fl_constant import FLContextKey
from nvflare.apis.fl_context import FLContext
from nvflare.app_common.aggregators.model_aggregator import ModelAggregator
from nvflare.client import FLModel


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TOOLS_ROOT = "/home/ubuntu/tools"
WRAPPER_NAME = "run_gwama.sh"
COMPARE_SCRIPT = "compare_gwama_models.py"
PLOT_SCRIPT = "plot_gwama_manhattan.py"
DASHBOARD_BUILDER = "build_fedx_dashboard.py"
DASHBOARD_HTML = "fedx_dashboard.html"

# The site driver currently implements REGENIE. Add another method here only
# when local_script_start_gwas.sh has a working branch and a tested converter.
METHODS = ("regenie",)
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


def _safe_component(value):
    """Make site-controlled identifiers safe for use in output filenames."""
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._")
    return cleaned or "unknown"


def _run_checked(cmd, label, env=None):
    """Run one central step and preserve its stdout/stderr in the job log."""
    print(f"\n{label}")
    print("Command: " + " ".join(map(str, cmd)))
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"{label} exited {result.returncode}")


class GWASMetaAggregator(ModelAggregator):
    """Collect GWAMA-format summary statistics and meta-analyse them."""

    def __init__(self, trait_type="binary", model="both", method="regenie",
                 tools_root=DEFAULT_TOOLS_ROOT, output_folder="server_results",
                 make_plots=True, make_dashboard=False):
        super().__init__()
        if trait_type not in TRAIT_TYPES:
            raise ValueError(f"unsupported trait_type: {trait_type}")
        if model not in ("fixed", "random", "both"):
            raise ValueError(f"unsupported meta-analysis model: {model}")
        if make_plots and model != "both":
            raise ValueError("matched Manhattan plots require model='both'")
        self.trait_type = trait_type
        self.meta_model = model
        self.method = method
        self.tools_root = tools_root
        self.output_folder = output_folder
        self.make_plots = make_plots
        self.make_dashboard = make_dashboard

        self.received_params_type = None
        self.passthrough_params = None
        self.output_dir = None
        self.accepted = []          # (dataset_id, site_name, path)
        self.rejected = []          # (site_name, reason)
        self.accepted_ids = set()
        self.genome_build = None

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

        site_method = meta.get("method")
        if site_method and site_method != self.method:
            return reject(f"ran method '{site_method}', "
                          f"server expects '{self.method}'")

        site_build = meta.get("build")
        if self.genome_build and site_build and site_build != self.genome_build:
            return reject(f"uses genome build '{site_build}', "
                          f"other sites use '{self.genome_build}'")

        if dataset_id in self.accepted_ids:
            return reject(f"duplicate dataset_id '{dataset_id}'")

        safe_id = _safe_component(dataset_id)
        safe_site = _safe_component(site)
        path = os.path.join(self.output_dir,
                            f"site{safe_id}_{safe_site}_gwama.txt")
        with open(path, "w") as f:
            f.write(content)

        qc = meta.get("qc_summary", "")
        if qc:
            with open(os.path.join(self.output_dir,
                                   f"site{safe_id}_{safe_site}_qc.txt"), "w") as f:
                f.write(qc)

        problem = self._check(path)
        if problem:
            return reject(problem)

        if self.genome_build is None and site_build:
            self.genome_build = site_build
        self.accepted.append((dataset_id, site, path))
        self.accepted_ids.add(dataset_id)
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
        files = [path for _, _, path in sorted(self.accepted)]
        cmd = ["bash", wrapper, mode, "--model", self.meta_model,
               output_root, *files]

        env = os.environ.copy()
        env["FEDGX_TOOLS_ROOT"] = self.tools_root

        _run_checked(cmd, f"Meta-analysis over {len(files)} sites", env=env)

        if self.meta_model == "both":
            meta_outputs = {
                "fixed": f"{output_root}.fixed.out",
                "random": f"{output_root}.random.out",
            }
            artifacts = {
                **meta_outputs,
                "comparison": f"{output_root}.comparison.tsv",
            }
        else:
            meta_outputs = {self.meta_model: f"{output_root}.out"}
            artifacts = dict(meta_outputs)

        missing = [path for path in artifacts.values()
                   if not os.path.isfile(path) or os.path.getsize(path) == 0]
        if missing:
            raise RuntimeError(f"central analysis did not produce: {missing}")

        helper_dir = os.path.dirname(wrapper)
        if self.make_plots:
            plot_script = os.path.join(helper_dir, PLOT_SCRIPT)
            if not os.path.isfile(plot_script):
                raise FileNotFoundError(f"plot helper not found: {plot_script}")
            plot_cmd = [
                sys.executable, plot_script,
                "--fixed", meta_outputs["fixed"],
                "--random", meta_outputs["random"],
                "--output-prefix", output_root,
            ]
            _run_checked(plot_cmd, "Creating matched Manhattan plots")
            artifacts.update({
                "fixed_manhattan": f"{output_root}.fixed.manhattan.png",
                "random_manhattan": f"{output_root}.random.manhattan.png",
            })

        if self.make_dashboard:
            builder = os.path.join(helper_dir, DASHBOARD_BUILDER)
            html_template = os.path.join(helper_dir, DASHBOARD_HTML)
            for required in (builder, html_template):
                if not os.path.isfile(required):
                    raise FileNotFoundError(f"dashboard helper not found: {required}")

            for model_name, meta_path in meta_outputs.items():
                dashboard_dir = os.path.join(
                    self.output_dir, "dashboard", model_name
                )
                dashboard_cmd = [
                    sys.executable, builder,
                    "--meta", meta_path,
                    "--model", model_name,
                    "--trait-type", self.trait_type,
                    "--method", self.method,
                    "--html-template", html_template,
                    "--out-dir", dashboard_dir,
                ]
                for site_path in files:
                    dashboard_cmd.extend(["--site", site_path])
                _run_checked(
                    dashboard_cmd,
                    f"Building {model_name}-effects interactive dashboard",
                )
                artifacts[f"{model_name}_dashboard"] = os.path.join(
                    dashboard_dir, "index.html"
                )

        final_missing = [path for path in artifacts.values()
                         if not os.path.isfile(path) or os.path.getsize(path) == 0]
        if final_missing:
            raise RuntimeError(f"post-processing did not produce: {final_missing}")

        return FLModel(
            params=self.passthrough_params,
            params_type=self.received_params_type,
            meta={
                "META_ANALYSIS_COMPLETED": True,
                "N_SITES": len(files),
                "N_REJECTED": len(self.rejected),
                "ARTIFACTS": artifacts,
            },
        )

    def reset_stats(self):
        self.accepted = []
        self.accepted_ids = set()
        self.rejected = []
        self.received_params_type = None
        self.passthrough_params = None
        self.genome_build = None
