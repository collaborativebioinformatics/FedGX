#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Server side code based off
https://github.com/collaborativebioinformatics/FedGen/blob/main/jobs/fed_gwas/job.py

Standardisation happens at the SITES, in gwas_cohort_qc_with_gwama.R, which
emits GWAMA-format summary statistics directly. The server therefore validates
what it receives rather than converting it, and never needs to know which GWAS
tool produced the numbers.
"""

import argparse
import os

from nvflare.app_opt.pt.recipes.fedavg import FedAvgRecipe
from nvflare.recipe import SimEnv, add_experiment_tracking, ProdEnv

from nvflare.app_common.aggregators.model_aggregator import ModelAggregator
from nvflare.client import FLModel

from nvflare.apis.fl_constant import FLContextKey
from nvflare.apis.fl_context import FLContext

from serverSide_metaAnalysis import runGWAMA

# Which client script runs which tool. No converter here: the R QC script at
# the site produces GWAMA format for every tool.
ADAPTERS = {
    "regenie": {"script": "client_regenie.sh"},
    "saige": {"script": "client_saige.sh"},
}

# Columns the R script writes, by trait type. The server checks for these and
# rejects a non-conforming site rather than trying to repair it.
REQUIRED_COLUMNS = {
    "binary": ["MARKERNAME", "EA", "NEA", "EAF", "N", "OR", "OR_95L", "OR_95U"],
    "quantitative": ["MARKERNAME", "EA", "NEA", "EAF", "N", "BETA", "SE"],
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

    def __init__(self, trait_type="binary", output_folder="server_results"):
        super().__init__()
        self.received_params_type = None
        self.trait_type = trait_type
        self.output_folder = output_folder

        # Set on the first accepted model
        self.gwama_input_file = None
        self.job_dir = None
        self.output_dir = None
        self.gwama_files = []
        self.rejected = []

    # -- helpers ------------------------------------------------------------

    def _ensure_output_dir(self):
        if self.output_dir is not None:
            return
        self.job_dir = _get_run_dir(self.fl_ctx)
        self.output_dir = os.path.join(self.job_dir, self.output_folder)
        os.makedirs(self.output_dir, exist_ok=True)
        self.gwama_input_file = os.path.join(self.output_dir, "gwama.in")
        print(f"Output directory: {self.output_dir}")

    def _validate(self, path, site_name):
        """
        Structural checks on a site file. Returns None if usable, otherwise a
        string describing why it is not.
        """
        required = REQUIRED_COLUMNS[self.trait_type]

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

        site_name = model.meta.get("site_name", "unknown_site")
        dataset_id = str(model.meta.get("dataset_id", "unknown_id"))

        if model.params.get("SUCCESS") is False:
            error_msg = model.meta.get("error_message", "Unknown error")
            print(f"ERROR: {site_name} reported failure: {error_msg}")
            self.rejected.append((site_name, f"client error: {error_msg}"))
            return

        results_file_content = model.meta.get("results_file", "")
        if not results_file_content:
            print(f"WARNING: no results_file content from {site_name}")
            self.rejected.append((site_name, "empty results_file"))
            return

        gwama_path = os.path.join(
            self.output_dir, f"site{dataset_id}_{site_name}_gwama.txt"
        )
        with open(gwama_path, "w") as f:
            f.write(results_file_content)
        print(f"Received {len(results_file_content)} bytes from {site_name} "
              f"(site{dataset_id}) -> {gwama_path}")

        # Optional QC report from the R script, stored alongside for the record
        qc_summary = model.meta.get("qc_summary", "")
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
        Write gwama.in in a deterministic order and run the meta-analysis.
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

        # Sorted so GWAMA's reference allele -- taken from the first file in the
        # list -- does not depend on the order the sites happened to respond in.
        with open(self.gwama_input_file, "w") as f:
            for _, path in sorted(self.gwama_files):
                f.write(f"{path}\n")

        print(f"\ngwama.in: {len(self.gwama_files)} sites")
        for dataset_id, path in sorted(self.gwama_files):
            print(f"  site{dataset_id}: {os.path.basename(path)}")

        runGWAMA(
            self.gwama_input_file,
            os.path.join(self.output_dir, "gwama"),
            n_expected=len(self.gwama_files),
        )

        aggregated_params = {
            "META_ANALYSIS_COMPLETED": True,
            "N_SITES": len(self.gwama_files),
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
                        choices=sorted(REQUIRED_COLUMNS),
                        help="Trait type; must match the --trait-type given to "
                             "the site R script (default: %(default)s)")

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
        aggregator=GWASMetaAggregator(trait_type=args.trait_type),
    )
    add_experiment_tracking(recipe, tracking_type="tensorboard")

    # Ship the method-specific GWAS wrapper and the shared R QC script, so every
    # site standardises with identical code rather than a local copy.
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