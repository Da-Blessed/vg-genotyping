from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = ROOT / "scripts" / "validate_inputs.py"
REFERENCE = ROOT / "tests" / "data" / "reference.fa"
PANEL = ROOT / "tests" / "data" / "panel.vcf"


class ValidatorTests(unittest.TestCase):
    def test_vg_170_autoindex_uses_tmp_dir_option(self) -> None:
        rule = (ROOT / "workflow/rules/longread.smk").read_text()
        self.assertIn("--tmp-dir", rule)
        self.assertNotIn("--temp-dir", rule)

    def run_validator(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(VALIDATOR), *arguments],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )

    def test_panel_only(self) -> None:
        result = self.run_validator(
            "--reference", str(REFERENCE), "--vcf", str(PANEL),
            "--panel-only", "--strict",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("2 variants", result.stdout)

    def test_short_reads(self) -> None:
        result = self.run_validator(
            "--reference", str(REFERENCE), "--vcf", str(PANEL),
            "--samples", str(ROOT / "tests/data/samples.shortread.tsv"),
            "--read-type", "short", "--strict",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("short-read target samples", result.stdout)

    def test_long_reads(self) -> None:
        result = self.run_validator(
            "--reference", str(REFERENCE), "--vcf", str(PANEL),
            "--samples", str(ROOT / "tests/data/samples.longread.tsv"),
            "--read-type", "long", "--strict",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("long-read target samples", result.stdout)

    def test_invalid_long_read_platform(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            sheet = Path(directory) / "samples.tsv"
            sheet.write_text(
                "sample\tread\tplatform\n"
                f"bad\t{ROOT / 'tests/data/long.fastq'}\tpacbio\n"
            )
            result = self.run_validator(
                "--reference", str(REFERENCE), "--vcf", str(PANEL),
                "--samples", str(sheet), "--read-type", "long", "--strict",
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("hifi or r10", result.stderr)


if __name__ == "__main__":
    unittest.main()
