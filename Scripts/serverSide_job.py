#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Server side code based off https://github.com/collaborativebioinformatics/FedGen/blob/main/jobs/fed_gwas/job.py
"""

import argparse
import os
from nvflare.app_opt.pt.recipes.fedavg import FedAvgRecipe
from nvflare.recipe import SimEnv, add_experiment_tracking, ProdEnv

from nvflare.app_common.aggregators.model_aggregator import ModelAggregator
from nvflare.client import FLModel

from nvflare.apis.fl_constant import FLContextKey
from nvflare.apis.fl_context import FLContext

from regenie2gwama import regenie2gwama
from saige2gwama import saige2gwama

from serverSide_metaAnalysis import runGWAMA

ADAPTERS = {
    "regenie": {
        "script": "client_regenie.sh",
        "raw_suffix": "regenie_step2_Phen1.regenie",
        "convert": lambda src, dst: regenie2gwama(src, dst, mode="or"),
    },
    "saige": {
        "script": "client_saige.sh",
        "raw_suffix": "saige_step2.txt",
        "convert": saige2gwama,
    },
}


def _get_run_dir(fl_ctx: FLContext):
    job_id = fl_ctx.get_job_id()
    workspace = fl_ctx.get_prop(FLContextKey.WORKSPACE_OBJECT)
    run_dir = workspace.get_run_dir(job_id)
    return run_dir


class GWASMetaAggregator(ModelAggregator):
    """
    Collects per-site GWAS summary statistics, converts each to GWAMA input
    format, and hands the assembled file list to the meta-analysis step.
    """

    def __init__(self, method="regenie", output_folder="server_results"):
        super().__init__()
        self.received_params_type = None
        self.method = method
        self.output_folder = output_folder

        # Set on the first accepted model
        self.gwama_input_file = None
        self.job_dir = None
        self.output_dir = None
        self.gwama_files = []

    def accept_model(self, model: FLModel):
        """
        Called once per client. Saves the raw tool output, converts it to
        GWAMA format, and records the converted path.
        """
        if self.received_params_type is None:
            self.received_params_type = model.params_type

        params = model.params
        adapter = ADAPTERS[self.method]

        # Create output directory if it doesn't exist
        if self.output_dir is None:
            self.job_dir = _get_run_dir(self.fl_ctx)
            self.output_dir = os.path.join(self.job_dir, self.output_folder)
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
                print(f"Created output directory: {self.output_dir}")

            self.gwama_input_file = os.path.join(self.output_dir, "gwama.in")
            print(f"GWAMA input file will be written at: {self.gwama_input_file}")

        # Check if this is an error response
        if params.get("SUCCESS") is False:
            site_name = model.meta.get("site_name", "unknown_site")
            error_msg = model.meta.get("error_message", "Unknown error")
            print(f"ERROR: Client {site_name} failed with error: {error_msg}")
            return

        # Extract metadata
        site_name = model.meta.get("site_name", "unknown_site")
        dataset_id = model.meta.get("dataset_id", "unknown_id")
        results_file_content = model.meta.get("results_file", "")

        if not results_file_content:
            print(f"WARNING: No results_file content received from {site_name}")
            return

        # Save the raw tool output
        output_filename = f"site{dataset_id}_{site_name}_{adapter['raw_suffix']}"
        output_path = os.path.join(self.output_dir, output_filename)
        with open(output_path, "w") as f:
            f.write(results_file_content)

        # Convert to GWAMA input format
        gwama_path = os.path.join(
            self.output_dir, f"site{dataset_id}_{site_name}_gwama.txt"
        )
        adapter["convert"](output_path, gwama_path)
        self.gwama_files.append((str(dataset_id), gwama_path))

        print(f"Saved results from {site_name} (site{dataset_id}) to: {output_path}")
        print(f"File size: {len(results_file_content)} bytes")
        print(f"Converted to GWAMA format: {gwama_path}")

    def aggregate_model(self) -> FLModel:
        """
        Write gwama.in in a deterministic order and run the meta-analysis.
        """
        if not self.gwama_files:
            raise RuntimeError("no site results converted; nothing to meta-analyse")

        # Sorted so GWAMA's reference allele (taken from the first file) does
        # not depend on the order sites happened to respond in.
        with open(self.gwama_input_file, "w") as f:
            for _, path in sorted(self.gwama_files):
                f.write(f"{path}\n")
        print(f"gwama.in: {len(self.gwama_files)} sites")

        # RUN GWAMA
        runGWAMA(
            self.gwama_input_file,
            os.path.join(self.output_dir, "gwama"),
            n_expected=len(self.gwama_files),
        )

        aggregated_params = {
            "META_ANALYSIS_COMPLETED": True,
            "N_SITES": len(self.gwama_files),
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
        self.received_params_type = None


def define_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_clients", type=int, default=3)
    parser.add_argument("--num_rounds", type=int, default=1)
    parser.add_argument("--env", type=str, default="sim", choices=["sim", "prod"],
                        help="Environment to run in: 'sim' for SimEnv or 'prod' for ProdEnv "
                             "(default: %(default)s)")
    parser.add_argument("--startup_kit", type=str, default=os.environ.get("FEDGX_STARTUP_KIT"),
                        help="Startup kit path for prod (or set FEDGX_STARTUP_KIT)")
    parser.add_argument("--username", default=os.environ.get("FEDGX_USERNAME"),
                        help="Admin username for prod (or set FEDGX_USERNAME)")
    parser.add_argument("--method", type=str, default="regenie",
                        choices=sorted(ADAPTERS),
                        help="GWAS tool the clients run (default: %(default)s)")

    args = parser.parse_args()
    if args.env == "prod" and not (args.startup_kit and args.username):
        parser.error("--env prod requires --startup_kit and --username "
                     "(or FEDGX_STARTUP_KIT / FEDGX_USERNAME)")
    return args


def main():
    args = define_parser()

    n_clients = args.n_clients
    num_rounds = args.num_rounds

    recipe = FedAvgRecipe(
        name="fed_gwas",
        min_clients=n_clients,
        num_rounds=num_rounds,
        train_script="client.py",
        train_args=f"--method {args.method}",
        aggregator=GWASMetaAggregator(method=args.method),
    )
    add_experiment_tracking(recipe, tracking_type="tensorboard")

    # Send the method-specific GWAS script to all clients
    recipe.job.to_clients(ADAPTERS[args.method]["script"])

    # Select environment based on command-line argument
    if args.env == "sim":
        print(f"Using Simulation Environment with {n_clients} clients")
        env = SimEnv(num_clients=n_clients)
    else:
        print("Using Production Environment")
        print(f"  Startup kit: {args.startup_kit}")
        print(f"  Username: {args.username}")
        env = ProdEnv(startup_kit_location=args.startup_kit, username=args.username,
                      login_timeout=300)

    run = recipe.execute(env)
    print()
    print("Job Status is:", run.get_status())
    print("Result can be found in :", run.get_result())
    print()


if __name__ == "__main__":
    main()