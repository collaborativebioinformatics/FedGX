import importlib.util
import shutil
import sys
import types
import unittest
import uuid
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1]
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


class FakeFLModel:
    def __init__(self, params=None, params_type=None, meta=None):
        self.params = params
        self.params_type = params_type
        self.meta = meta or {}


def stub_module(name, **attributes):
    module = types.ModuleType(name)
    module.__dict__.update(attributes)
    sys.modules[name] = module
    return module


for package in (
    "nvflare",
    "nvflare.apis",
    "nvflare.app_common",
    "nvflare.app_common.aggregators",
    "nvflare.app_opt",
    "nvflare.app_opt.pt",
    "nvflare.app_opt.pt.recipes",
):
    stub_module(package)


stub_module("nvflare.apis.fl_constant", FLContextKey=types.SimpleNamespace(WORKSPACE_OBJECT="workspace"))
stub_module("nvflare.apis.fl_context", FLContext=object)
stub_module(
    "nvflare.app_common.aggregators.model_aggregator",
    ModelAggregator=type("ModelAggregator", (), {"__init__": lambda self: None}),
)
stub_module("nvflare.app_opt.pt.recipes.fedavg", FedAvgRecipe=object)
stub_module("nvflare.client", FLModel=FakeFLModel)
stub_module(
    "nvflare.recipe",
    ProdEnv=object,
    SimEnv=object,
    add_experiment_tracking=lambda *args, **kwargs: None,
)
stub_module("model", DummyModel=object)

SPEC = importlib.util.spec_from_file_location(
    "fedgx_meta_aggregator_test_module", SCRIPTS / "fedgx_meta_aggregator.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class ServerPipelineTests(unittest.TestCase):
    def setUp(self):
        self.root = Path.cwd() / ".test-output" / uuid.uuid4().hex
        self.root.mkdir(parents=True)

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        try:
            self.root.parent.rmdir()
        except OSError:
            pass

    def test_both_models_plots_and_optional_dashboards_are_connected(self):
        wrapper = self.root / MODULE.WRAPPER_NAME
        wrapper.write_text("#!/usr/bin/env bash\n", encoding="utf-8")
        for name in (MODULE.PLOT_SCRIPT, MODULE.DASHBOARD_BUILDER, MODULE.DASHBOARD_HTML):
            (self.root / name).write_text("test helper\n", encoding="utf-8")

        site1 = self.root / "site1.txt"
        site2 = self.root / "site2.txt"
        site1.write_text("MARKERNAME\n1:100:A:C\n", encoding="utf-8")
        site2.write_text("MARKERNAME\n1:100:A:C\n", encoding="utf-8")

        aggregator = MODULE.GWASMetaAggregator(
            trait_type="binary",
            model="both",
            method="regenie",
            make_plots=True,
            make_dashboard=True,
        )
        aggregator.output_dir = str(self.root)
        aggregator.accepted = [
            ("1", "site1", str(site1)),
            ("2", "site2", str(site2)),
        ]
        aggregator.passthrough_params = {"placeholder": 0}
        aggregator.received_params_type = "FULL"

        commands = []
        original_find = MODULE._find_wrapper
        original_run = MODULE._run_checked

        def fake_run(cmd, label, env=None):
            commands.append((cmd, label, env))
            output_root = self.root / "gwama"
            if cmd[0] == "bash":
                (Path(str(output_root) + ".fixed.out")).write_text("fixed\n", encoding="utf-8")
                (Path(str(output_root) + ".random.out")).write_text("random\n", encoding="utf-8")
                (Path(str(output_root) + ".comparison.tsv")).write_text("comparison\n", encoding="utf-8")
            elif Path(cmd[1]).name == MODULE.PLOT_SCRIPT:
                Path(str(output_root) + ".fixed.manhattan.png").write_bytes(b"fixed png")
                Path(str(output_root) + ".random.manhattan.png").write_bytes(b"random png")
            elif Path(cmd[1]).name == MODULE.DASHBOARD_BUILDER:
                dashboard_dir = Path(cmd[cmd.index("--out-dir") + 1])
                dashboard_dir.mkdir(parents=True, exist_ok=True)
                (dashboard_dir / "index.html").write_text("dashboard\n", encoding="utf-8")

        try:
            MODULE._find_wrapper = lambda: str(wrapper)
            MODULE._run_checked = fake_run
            result = aggregator.aggregate_model()
        finally:
            MODULE._find_wrapper = original_find
            MODULE._run_checked = original_run

        self.assertTrue(result.meta["META_ANALYSIS_COMPLETED"])
        self.assertEqual(result.meta["N_SITES"], 2)
        self.assertEqual(len(commands), 4)
        self.assertEqual(commands[0][0][0:5], ["bash", str(wrapper), "or", "--model", "both"])
        self.assertIn("fixed_manhattan", result.meta["ARTIFACTS"])
        self.assertIn("random_manhattan", result.meta["ARTIFACTS"])
        self.assertIn("fixed_dashboard", result.meta["ARTIFACTS"])
        self.assertIn("random_dashboard", result.meta["ARTIFACTS"])

    def test_filename_components_cannot_escape_output_directory(self):
        self.assertEqual(MODULE._safe_component("../../site one"), "site_one")

    def test_client_failure_metadata_is_rejected(self):
        aggregator = MODULE.GWASMetaAggregator()
        aggregator.output_dir = str(self.root)
        aggregator.accept_model(
            FakeFLModel(
                params={"placeholder": 0},
                params_type="FULL",
                meta={
                    "success": False,
                    "site_name": "site-1",
                    "dataset_id": "site1",
                    "error_message": "local GWAS failed",
                },
            )
        )
        self.assertEqual(aggregator.rejected, [("site-1", "local GWAS failed")])
        self.assertEqual(aggregator.accepted, [])


if __name__ == "__main__":
    unittest.main()
