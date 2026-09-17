import csv
import importlib.util
import shutil
import unittest
import uuid
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "compare_gwama_models.py"
SPEC = importlib.util.spec_from_file_location("compare_gwama_models", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class CompareGwamaModelsTests(unittest.TestCase):
    def setUp(self):
        self.root = Path.cwd() / ".test-output" / uuid.uuid4().hex
        self.root.mkdir(parents=True)

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        try:
            self.root.parent.rmdir()
        except OSError:
            pass

    def test_compares_binary_outputs_and_retains_unmatched_markers(self):
        header = (
            "rs_number reference_allele other_allele eaf OR OR_se OR_95L "
            "OR_95U z p-value q_statistic q_p-value i2 n_studies\n"
        )
        fixed = self.root / "fixed.out"
        random = self.root / "random.out"
        output = self.root / "comparison.tsv"
        fixed.write_text(
            header
            + "rs1 A C 0.2 1.10 0.02 1.06 1.14 4.0 0.001 3.2 0.20 30 3\n"
            + "rs2 G T 0.3 0.90 0.03 0.84 0.96 -3.0 0.003 1.0 0.60 0 3\n",
            encoding="utf-8",
        )
        random.write_text(
            header
            + "rs1 A C 0.2 1.08 0.05 0.98 1.18 1.5 0.10 3.2 0.20 30 3\n"
            + "rs3 C T 0.1 1.20 0.08 1.04 1.36 2.1 0.04 2.0 0.30 10 2\n",
            encoding="utf-8",
        )

        count = MODULE.compare(fixed, random, output)

        self.assertEqual(count, 3)
        with output.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual([row["MARKERNAME"] for row in rows], ["rs1", "rs2", "rs3"])
        self.assertEqual(rows[0]["STATUS"], "both")
        self.assertEqual(rows[0]["FFX_EFFECT"], "1.10")
        self.assertEqual(rows[0]["RFX_EFFECT"], "1.08")
        self.assertEqual(rows[0]["EFFECT_ABS_DIFF"], "0.02")
        self.assertEqual(rows[0]["Q_STATISTIC"], "3.2")
        self.assertEqual(rows[1]["STATUS"], "fixed_only")
        self.assertEqual(rows[2]["STATUS"], "random_only")

    def test_compares_quantitative_outputs(self):
        header = "rs_number reference_allele other_allele beta se p-value n_studies\n"
        fixed = self.root / "fixed.out"
        random = self.root / "random.out"
        output = self.root / "comparison.tsv"
        fixed.write_text(header + "rs1 A C 0.10 0.02 1e-6 3\n", encoding="utf-8")
        random.write_text(header + "rs1 A C 0.08 0.05 0.11 3\n", encoding="utf-8")

        MODULE.compare(fixed, random, output)

        with output.open(encoding="utf-8", newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["FFX_EFFECT"], "0.10")
        self.assertEqual(row["RFX_SE"], "0.05")

    def test_rejects_duplicate_markers(self):
        path = self.root / "duplicate.out"
        path.write_text(
            "rs_number OR OR_se p-value\n"
            "rs1 1.1 0.1 0.2\n"
            "rs1 1.2 0.1 0.1\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "duplicate marker rs1"):
            MODULE.read_gwama(path)


if __name__ == "__main__":
    unittest.main()
