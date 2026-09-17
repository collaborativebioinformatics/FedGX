import importlib.util
import shutil
import unittest
import uuid
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "plot_gwama_manhattan.py"
SPEC = importlib.util.spec_from_file_location("plot_gwama_manhattan", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class PlotGwamaManhattanTests(unittest.TestCase):
    def setUp(self):
        self.root = Path.cwd() / ".test-output" / uuid.uuid4().hex
        self.root.mkdir(parents=True)

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        try:
            self.root.parent.rmdir()
        except OSError:
            pass

    def write_result(self, path, rows):
        path.write_text(
            "rs_number p-value _-log10_p-value\n" + "\n".join(rows) + "\n",
            encoding="utf-8",
        )

    def test_marker_parser(self):
        self.assertEqual(MODULE.parse_marker("chr01:12345:A:C"), ("1", 12345))
        self.assertEqual(MODULE.parse_marker("X:50:G:T"), ("X", 50))
        with self.assertRaises(ValueError):
            MODULE.parse_marker("rs123")

    def test_creates_two_matched_manhattan_plots(self):
        fixed = self.root / "meta.fixed.out"
        random = self.root / "meta.random.out"
        prefix = self.root / "meta"
        self.write_result(
            fixed,
            [
                "1:100:A:C 1e-9 9",
                "1:200:G:T 0.01 2",
                "2:50:C:T 1e-6 6",
            ],
        )
        self.write_result(
            random,
            [
                "1:100:A:C 1e-4 4",
                "1:200:G:T 0.02 1.69897",
                "2:50:C:T 1e-3 3",
            ],
        )

        fixed_plot, random_plot, fixed_skipped, random_skipped = MODULE.create_plots(
            fixed, random, prefix, dpi=80
        )

        self.assertEqual((fixed_skipped, random_skipped), (0, 0))
        self.assertGreater(fixed_plot.stat().st_size, 1000)
        self.assertGreater(random_plot.stat().st_size, 1000)

    def test_requires_coordinate_marker_names(self):
        result = self.root / "meta.out"
        self.write_result(result, ["rs123 0.01 2"])
        with self.assertRaisesRegex(ValueError, "CHR:POS"):
            MODULE.read_gwama(result)


if __name__ == "__main__":
    unittest.main()
