#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Server side code based off
https://github.com/collaborativebioinformatics/FedGen/blob/main/jobs/fed_gwas/job.py

Standardisation happens at the SITES, in gwas_cohort_qc_with_gwama.R, which
emits GWAMA-format summary statistics directly. The server therefore validates
what it receives rather than converting it, and never needs to know which GWAS
tool produced the numbers.

The meta-analysis itself is delegated to run_gwama.sh, so there is exactly one
GWAMA invocation in the codebase: the same script can be re-run by hand against
the collected files without repeating the GWAS.
"""

import argparse
import os
import subprocess

from nvflare.app_opt.pt.recipes.fedavg import FedAvgRecipe
from nvflare.recipe import SimEnv, add_experiment_tracking, ProdEnv

from nvflare.app_common.aggregators.model_aggregator import ModelAggregator
from nvflare.client import FLModel

from nvflare.apis.fl_constant import FLContextKey
from nvflare.apis.fl_context import FLContext

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GWAMA_WRAPPER = os.environ.get(
    "FEDGX_GWAMA_WRAPPER", os.path.join(SCRIPT_DIR, "run_gwama.sh")
)

# Every method uses the same site driver, which dispatches internally on
# PROGRAM. Listing them separately keeps --method validated and gives a place
# to hang per-method differences later.
ADAPTERS = {
    "regenie": {"script": "local_script_start_gwas.sh"},
    "saige": {"script": "local_script_start_gwas.sh"},
}

# Columns the R script writes, by trait type, and the mode run_gwama.sh needs.
TRAIT_TYPES = {
    "binary": {
        "gwama_mode": "or",
        "columns": ["MARKERNAME", "EA", "NEA", "EAF", "N",
                    "OR", "OR_95L", "OR_95U"],
    },
    "quantitative": {
        "gwama_mode": "qt",
        "columns": ["MARKERNAME", "EA", "NEA", "EAF", "N", "BETA", "SE"],
    },
}


def _get_run_dir(fl_ctx: FLContext):
    job_id = fl_ctx.get_job_id()
    workspace = fl_ctx.get_prop(FLContextKey.WORKSPACE_OBJECT)
    run_dir = workspace.get_run_dir(job_id)
    return run_dir


class GWASMetaAggregator(ModelAggregator):
    """
    Collects GWAMA-format summary statistics from the sites, validates each
    file, and runs the meta-analysis over the ones that pass.
    """

    def __init__(self, trait_type="binary", model="fixed",
                 output_folder="server_results"):
        super().__init__()
        self.received_params_type = None
        self.trait_type = trait_type
        self.model = model
        self.output_folder = output_folder

        # Set on the first accepted model
        self.job_dir = None
        self.output_dir = None
        self.gwama_files = []
        self.rejected = []

    # -- helpers -------------------------------------------------------------

    def _ensure_output_dir(self):
        if self.output_dir is not None:
            return
        self.job_dir = _get_run_dir(self.fl_ctx)
        self.output_dir = os.path.join(self.job_dir, self.output_folder)
        os.makedirs(self.output_dir, exist_ok=True)
        print(f"Output directory: {self.output_dir}")

    def _validate(self, path, site_name):
        """
        Structural checks on a site file. Returns None if usable, otherwise a
        string describing why it is not.
        """
        required = TRAIT_TYPES[self.trait_type]["columns"]

        try:
            with open(path) as f:
                header_line = f.readline()
                if not header_line:
                    return "file is empty"
                header = header_line.rstrip("\n").split("\t")
                n_rows = sum(1 for _ in f)
        except OSError as e:
            return f"could not read file: {e}"

        missing = [c for c in required if c not in header]
        if missing:
            return (f"missing column(s) {missing} for trait type "
                    f"'{self.trait_type}'; header was {header}")

        if n_rows == 0:
            return "no data rows"

        print(f"  {site_name}: {n_rows} variants, columns OK")
        return None

    # -- aggregator interface ------------------------------------------------

    def accept_model(self, model: FLModel):
        """
        Called once per client. The payload is already GWAMA-format text from
        the site's R QC step, so the server only stores and checks it.
        """
        if self.received_params_type is None:
            self.received_params_type = model.params_type

        self._ensure_output_dir()

        meta = model.meta or {}
        site_name = meta.get("site_name", "unknown_site")
        dataset_id = str(meta.get("dataset_id", "unknown_id"))

        if model.params.get("SUCCESS") is False:
            error_msg = meta.get("error_message", "Unknown error")
            print(f"ERROR: {site_name} reported failure: {error_msg}")
            self.rejected.append((site_name, f"client error: {error_msg}"))
            return

        results_file_content = meta.get("results_file", "")
        if not results_file_content:
            print(f"WARNING: no results_file content from {site_name}")
            self.rejected.append((site_name, "empty results_file"))
            return

        # A site running a different trait type would produce the wrong columns
        # anyway, but saying so explicitly gives a clearer error.
        site_trait = meta.get("trait_type")
        if site_trait and site_trait != self.trait_type:
            print(f"ERROR: {site_name} ran trait_type '{site_trait}', "
                  f"server expects '{self.trait_type}'")
            self.rejected.append(
                (site_name, f"trait_type mismatch: {site_trait}")
            )
            return

        gwama_path = os.path.join(
            self.output_dir, f"site{dataset_id}_{site_name}_gwama.txt"
        )
        with open(gwama_path, "w") as f:
            f.write(results_file_content)
        print(f"Received {len(results_file_content)} bytes from {site_name} "
              f"(site{dataset_id}) -> {gwama_path}")

        qc_summary = meta.get("qc_summary", "")
        if qc_summary:
            qc_path = os.path.join(
                self.output_dir, f"site{dataset_id}_{site_name}_qc_summary.txt"
            )
            with open(qc_path, "w") as f:
                f.write(qc_summary)
            print(f"  QC summary saved to {qc_path}")

        problem = self._validate(gwama_path, site_name)
        if problem:
            print(f"REJECTED {site_name}: {problem}")
            self.rejected.append((site_name, problem))
            return

        self.gwama_files.append((dataset_id, gwama_path))

    def aggregate_model(self) -> FLModel:
        """
        Hand the collected files to run_gwama.sh, which sorts the filelist,
        validates the headers, runs GWAMA and checks n_studies afterwards.
        """
        if self.rejected:
            print("\nSites excluded from the meta-analysis:")
            for site_name, reason in self.rejected:
                print(f"  {site_name}: {reason}")

        if len(self.gwama_files) < 2:
            raise RuntimeError(
                f"only {len(self.gwama_files)} usable site file(s); "
                "a meta-analysis needs at least 2"
            )

        if not os.path.isfile(GWAMA_WRAPPER):
            raise FileNotFoundError(
                f"GWAMA wrapper not found: {GWAMA_WRAPPER} "
                "(set FEDGX_GWAMA_WRAPPER)"
            )

        mode = TRAIT_TYPES[self.trait_type]["gwama_mode"]
        output_root = os.path.join(self.output_dir, "gwama")

        # Sorted here as well as inside the wrapper, so the command that gets
        # logged is the command that can be replayed verbatim.
        files = [path for _, path in sorted(self.gwama_files)]

        cmd = ["bash", GWAMA_WRAPPER, mode, "--model", self.model,
               output_root] + files

        print(f"\nMeta-analysis over {len(files)} sites")
        print("Command: " + " ".join(cmd))

        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(result.stderr)

        if result.returncode != 0:
            raise RuntimeError(
                f"run_gwama.sh exited {result.returncode}; see output above"
            )

        aggregated_params = {
            "META_ANALYSIS_COMPLETED": True,
            "N_SITES": len(files),
            "N_REJECTED": len(self.rejected),
        }

        return FLModel(
            params=aggregated_params,
            params_type=self.received_params_type,
        )

    def reset_stats(self):
        """
        Clear state between FL rounds.
        """
        self.gwama_files = []
        self.rejected = []
        self.received_params_type = None


def define_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_clients", type=int, default=3)
    parser.add_argument("--num_rounds", type=int, default=1)
    parser.add_argument("--env", type=str, default="sim", choices=["sim", "prod"],
                        help="Environment to run in: 'sim' for SimEnv or 'prod' "
                             "for ProdEnv (default: %(default)s)")
    parser.add_argument("--startup_kit", type=str,
                        default=os.environ.get("FEDGX_STARTUP_KIT"),
                        help="Startup kit path for prod (or set FEDGX_STARTUP_KIT)")
    parser.add_argument("--username", default=os.environ.get("FEDGX_USERNAME"),
                        help="Admin username for prod (or set FEDGX_USERNAME)")
    parser.add_argument("--method", type=str, default="regenie",
                        choices=sorted(ADAPTERS),
                        help="GWAS tool the clients run (default: %(default)s)")
    parser.add_argument("--trait_type", type=str, default="binary",
                        choices=sorted(TRAIT_TYPES),
                        help="Trait type; the clients pass this to the R QC "
                             "script (default: %(default)s)")
    parser.add_argument("--model", type=str, default="fixed",
                        choices=["fixed", "random", "both"],
                        help="Meta-analysis model. With few cohorts, random "
                             "effects are unstable (default: %(default)s)")

    args = parser.parse_args()
    if args.env == "prod" and not (args.startup_kit and args.username):
        parser.error("--env prod requires --startup_kit and --username "
                     "(or FEDGX_STARTUP_KIT / FEDGX_USERNAME)")
    return args


def main():
    args = define_parser()

    recipe = FedAvgRecipe(
        name="fed_gwas",
        min_clients=args.n_clients,
        num_rounds=args.num_rounds,
        train_script="client.py",
        train_args=f"--method {args.method} --trait_type {args.trait_type}",
        aggregator=GWASMetaAggregator(trait_type=args.trait_type,
                                      model=args.model),
    )
    add_experiment_tracking(recipe, tracking_type="tensorboard")

    # Ship the site driver and the shared R QC script, so every site
    # standardises with identical code rather than a local copy.
    recipe.job.to_clients(ADAPTERS[args.method]["script"])
    recipe.job.to_clients("gwas_cohort_qc_with_gwama.R")

    if args.env == "sim":
        print(f"Using Simulation Environment with {args.n_clients} clients")
        env = SimEnv(num_clients=args.n_clients)
    else:
        print("Using Production Environment")
        print(f"  Startup kit: {args.startup_kit}")
        print(f"  Username: {args.username}")
        env = ProdEnv(startup_kit_location=args.startup_kit,
                      username=args.username, login_timeout=300)

    run = recipe.execute(env)
    print()
    print("Job Status is:", run.get_status())
    print("Result can be found in :", run.get_result())
    print()


if __name__ == "__main__":
    main()