import subprocess
import sys
from pathlib import Path

import pytest

from hydrogen_s3.reference_spectroscopy import conventional_rydberg_baseline, load_reference_lines


def test_reference_lines_load_with_medium_labels() -> None:
    rows = load_reference_lines()
    assert rows
    assert {row.medium for row in rows} == {"air", "vacuum"}
    assert all(row.source_url for row in rows)


def test_vacuum_residuals_only_when_reference_medium_matches() -> None:
    rows = conventional_rydberg_baseline()
    assert rows
    for row in rows:
        if row["reference_medium"] == "vacuum":
            assert row["comparison_valid"] is True
            assert row["residual_nm"] is not None
            assert row["relative_error_ppm"] is not None
        else:
            assert row["comparison_valid"] is False
            assert row["residual_nm"] is None
            assert row["relative_error_ppm"] is None


def test_missing_reference_columns_are_rejected(tmp_path: Path) -> None:
    bad_csv = tmp_path / "bad.csv"
    bad_csv.write_text("series,transition\nLyman,2->1\n", encoding="utf-8")
    with pytest.raises(ValueError):
        load_reference_lines(bad_csv)


def test_report_generation() -> None:
    subprocess.run([sys.executable, "scripts/generate_report.py"], check=True)
    report = Path("build/reports/hydrogen_s3_report.md")
    assert report.exists()
    text = report.read_text(encoding="utf-8")
    assert "Branching Transform Diagnostics" in text
    assert "Conventional Rydberg Baseline" in text
