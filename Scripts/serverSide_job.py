#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Build and execute the FedGX NVFLARE job."""

import argparse
import os

from nvflare.app_opt.pt.recipes.fedavg import FedAvgRecipe
from nvflare.recipe import ProdEnv, SimEnv, add_experiment_tracking

from fedgx_meta_aggregator import (
    COMPARE_SCRIPT,
    DASHBOARD_BUILDER,
    DASHBOARD_HTML,
    DEFAULT_TOOLS_ROOT,
    GWASMetaAggregator,
    METHODS,
    PLOT_SCRIPT,
    QC_SCRIPT,
    SITE_SCRIPT,
    TRAIT_TYPES,
    WRAPPER_NAME,
)
from model import DummyModel


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_SCRIPT = "model.py"
AGGREGATOR_SCRIPT = "fedgx_meta_aggregator.py"


def script_path(name):
    """Return a launch-directory-independent source path for job packaging."""
    return os.path.join(SCRIPT_DIR, name)


def define_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Run one federated GWAS round, central GWAMA meta-analysis, "
            "matched FFX/RFX plots, and an optional interactive dashboard."
        )
    )
    parser.add_argument("--n_clients", type=int, default=3)
    parser.add_argument("--num_rounds", type=int, default=1)
    parser.add_argument(
        "--env", default="sim", choices=["sim", "prod"],
        help="default: %(default)s",
    )
    parser.add_argument(
        "--method", default="regenie", choices=METHODS,
        help="GWAS tool the clients run (default: %(default)s)",
    )
    parser.add_argument(
        "--trait_type", default="binary", choices=sorted(TRAIT_TYPES),
        help="passed to the site R script (default: %(default)s)",
    )
    parser.add_argument(
        "--model", default="both", choices=["fixed", "random", "both"],
        help="meta-analysis model (default: %(default)s)",
    )
    parser.add_argument(
        "--plots", action=argparse.BooleanOptionalAction, default=True,
        help=(
            "create matched FFX/RFX Manhattan PNGs; requires --model both "
            "(default: %(default)s)"
        ),
    )
    parser.add_argument(
        "--dashboard", action="store_true",
        help="also build optional fixed/random interactive dashboards",
    )
    parser.add_argument(
        "--tools_root",
        default=os.environ.get("FEDGX_TOOLS_ROOT", DEFAULT_TOOLS_ROOT),
        help="where GWAMA lives on the server (default: %(default)s)",
    )
    parser.add_argument(
        "--startup_kit", default=os.environ.get("FEDGX_STARTUP_KIT"),
        help="prod only (or set FEDGX_STARTUP_KIT)",
    )
    parser.add_argument(
        "--username", default=os.environ.get("FEDGX_USERNAME"),
        help="prod only (or set FEDGX_USERNAME)",
    )

    args = parser.parse_args()
    if args.env == "prod" and not (args.startup_kit and args.username):
        parser.error(
            "--env prod requires --startup_kit and --username "
            "(or FEDGX_STARTUP_KIT / FEDGX_USERNAME)"
        )
    if args.plots and args.model != "both":
        parser.error(
            "--plots requires --model both; use --no-plots for a single model"
        )
    if args.num_rounds != 1:
        parser.error(
            "FedGX performs one complete GWAS/meta-analysis round; "
            "--num_rounds must be 1"
        )
    return args


def build_recipe(args):
    """Construct the job and bundle every runtime dependency explicitly."""
    recipe = FedAvgRecipe(
        name="fed_gwas",
        model=DummyModel(),
        min_clients=args.n_clients,
        num_rounds=args.num_rounds,
        train_script=script_path("client.py"),
        train_args=f"--method {args.method} --trait_type {args.trait_type}",
        aggregator=GWASMetaAggregator(
            trait_type=args.trait_type,
            model=args.model,
            method=args.method,
            tools_root=args.tools_root,
            make_plots=args.plots,
            make_dashboard=args.dashboard,
        ),
    )
    add_experiment_tracking(recipe, tracking_type="tensorboard")

    # Client runtime: site GWAS driver and identical QC/conversion code.
    recipe.job.to_clients(script_path(SITE_SCRIPT))
    recipe.job.to_clients(script_path(QC_SCRIPT))

    # Server runtime: importable custom classes plus central analysis helpers.
    recipe.job.to_server(script_path(AGGREGATOR_SCRIPT))
    recipe.job.to_server(script_path(MODEL_SCRIPT))
    recipe.job.to_server(script_path(WRAPPER_NAME))
    recipe.job.to_server(script_path(COMPARE_SCRIPT))
    if args.plots:
        recipe.job.to_server(script_path(PLOT_SCRIPT))
    if args.dashboard:
        recipe.job.to_server(script_path(DASHBOARD_BUILDER))
        recipe.job.to_server(script_path(DASHBOARD_HTML))

    return recipe


def main():
    args = define_parser()
    recipe = build_recipe(args)

    if args.env == "sim":
        print(f"Simulation environment, {args.n_clients} clients")
        environment = SimEnv(num_clients=args.n_clients)
    else:
        print(f"Production environment, startup kit {args.startup_kit}")
        environment = ProdEnv(
            startup_kit_location=args.startup_kit,
            username=args.username,
            login_timeout=300,
        )

    run = recipe.execute(environment)
    print("\nJob status:", run.get_status())
    print("Results in:", run.get_result())


if __name__ == "__main__":
    main()
