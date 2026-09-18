import json
import math
import shutil
import subprocess
import sys
import unittest
import uuid
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1]
BUILDER = SCRIPTS / "build_fedx_dashboard.py"
HTML_TEMPLATE = SCRIPTS / "fedx_dashboard.html"


SITE_HEADER = "MARKERNAME\tCHR\tPOS\tEA\tNEA\tEAF\tN\tOR\tOR_95L\tOR_95U\n"
META_HEADER = (
    "rs_number reference_allele other_allele eaf OR OR_se OR_95L OR_95U "
    "z p-value _-log10_p-value q_statistic q_p-value i2 n_studies "
    "n_samples effects\n"
)


def read_payload(path):
    text = path.read_text(encoding="utf-8")
    prefix = "window.FEDX_DATA = "
    return json.loads(text[text.index(prefix) + len(prefix):].rstrip(";\n"))


class DashboardBuilderTests(unittest.TestCase):
    def setUp(self):
        self.root = Path.cwd() / ".test-output" / uuid.uuid4().hex
        self.root.mkdir(parents=True)

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        try:
            self.root.parent.rmdir()
        except OSError:
            pass

    def test_binary_gwama_inputs_are_normalized_to_log_effects(self):
        root = self.root
        site1 = root / "site1.GWAMA.txt"
        site2 = root / "site2.GWAMA.txt"
        meta = root / "gwama.fixed.out"
        output = root / "dashboard"

        site1.write_text(
            SITE_HEADER
            + "1:100:A:C\t1\t100\tA\tC\t0.20\t1000\t1.221403\t1.002167\t1.488609\n"
            + "2:200:G:T\t2\t200\tG\tT\t0.30\t1000\t0.904837\t0.743787\t1.100759\n",
            encoding="utf-8",
        )
        site2.write_text(
            SITE_HEADER
            + "1:100:A:C\t1\t100\tA\tC\t0.22\t1200\t1.161834\t0.990050\t1.363425\n"
            + "2:200:G:T\t2\t200\tG\tT\t0.31\t1200\t0.951229\t0.810584\t1.116278\n",
            encoding="utf-8",
        )
        meta.write_text(
            META_HEADER
            + "1:100:A:C A C 0.21 1.197217 0.05 1.085456 1.320488 3.6 0.00032 3.49 0.4 0.81 0 2 2200 ++\n"
            + "2:200:G:T G T 0.30 0.927743 0.05 0.841235 1.023146 -1.5 0.13 0.89 0.2 0.90 0 2 2200 --\n",
            encoding="utf-8",
        )

        completed = subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--site", str(site1),
                "--site", str(site2),
                "--meta", str(meta),
                "--model", "fixed",
                "--trait-type", "binary",
                "--method", "regenie",
                "--html-template", str(HTML_TEMPLATE),
                "--out-dir", str(output),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

        payload = read_payload(output / "fedx_data.js")
        self.assertEqual(payload["analysis"]["model"], "fixed")
        self.assertEqual(payload["analysis"]["traitType"], "binary")
        self.assertEqual(payload["analysis"]["method"], "REGENIE")
        self.assertEqual(len(payload["snps"]), 2)

        first = next(row for row in payload["snps"] if row["marker"] == "1:100:A:C")
        self.assertAlmostEqual(first["meta"]["beta"], math.log(1.197217), places=5)
        self.assertAlmostEqual(first["meta"]["z"], 3.6, places=6)
        self.assertAlmostEqual(first["meta"]["p"], 0.00032, places=8)
        self.assertAlmostEqual(first["site"]["site1"]["beta"], 0.2, places=5)

        html = (output / "index.html").read_text(encoding="utf-8")
        self.assertTrue(html.lower().startswith("<!doctype html>"))
        self.assertIn("dataModeBanner", html)

    def test_quantitative_random_effect_dashboard(self):
        site = self.root / "site1.GWAMA.txt"
        meta = self.root / "gwama.random.out"
        output = self.root / "random-dashboard"
        site.write_text(
            "MARKERNAME\tCHR\tPOS\tEA\tNEA\tEAF\tN\tBETA\tSE\n"
            "3:300:A:G\t3\t300\tA\tG\t0.25\t900\t0.12\t0.04\n"
            "4:400:C:T\t4\t400\tC\tT\t0.30\t900\t0.01\t0.05\n",
            encoding="utf-8",
        )
        meta.write_text(
            "rs_number beta beta_se z p-value n_samples\n"
            "3:300:A:G 0.10 0.05 2.0 0.0455 900\n"
            "4:400:C:T 0.01 0.05 0.2 0.84 900\n",
            encoding="utf-8",
        )

        completed = subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--site", str(site),
                "--meta", str(meta),
                "--model", "random",
                "--trait-type", "quantitative",
                "--max-variants", "1",
                "--out-dir", str(output),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = read_payload(output / "fedx_data.js")
        self.assertEqual(payload["analysis"]["model"], "random")
        self.assertAlmostEqual(payload["snps"][0]["meta"]["beta"], 0.10)
        self.assertIn("random-effects", payload["meta"]["sub"])
        self.assertEqual(payload["sampling"]["totalVariants"], 2)
        self.assertEqual(payload["sampling"]["displayedVariants"], 1)
        self.assertTrue(payload["sampling"]["sampled"])


if __name__ == "__main__":
    unittest.main()
