"""Conventional Rydberg baseline against auditable reference lines.

The comparison uses a vacuum shell-law calculation. Rows with non-vacuum
reference media keep their raw NIST wavelength and medium, but residual fields
are suppressed rather than mixing air and vacuum conventions.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from importlib.resources import files
from importlib.resources.abc import Traversable
from pathlib import Path

from hydrogen_s3.spectrum import vacuum_transition_wavelength_nm

REFERENCE_DATA_RESOURCE = files("hydrogen_s3").joinpath("data", "hydrogen_reference_lines.csv")

REFERENCE_COLUMNS = {
    "series",
    "transition",
    "n_i",
    "n_f",
    "wavelength_nm",
    "medium",
    "source",
    "source_url",
    "source_access_date",
    "source_table_or_query",
    "notes",
}


@dataclass(frozen=True)
class ReferenceLine:
    """Raw reference line metadata from the repository CSV fixture."""

    series: str
    transition: str
    n_i: int
    n_f: int
    wavelength_nm: float
    medium: str
    source: str
    source_url: str
    source_access_date: str
    source_table_or_query: str
    notes: str


def load_reference_lines(path: Path | Traversable | None = None) -> list[ReferenceLine]:
    """Load auditable reference rows from CSV."""

    data_path = REFERENCE_DATA_RESOURCE if path is None else path
    with data_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = REFERENCE_COLUMNS - fields
        if missing:
            raise ValueError(f"Missing required reference columns: {sorted(missing)}")
        rows = []
        for row in reader:
            rows.append(
                ReferenceLine(
                    series=row["series"],
                    transition=row["transition"],
                    n_i=int(row["n_i"]),
                    n_f=int(row["n_f"]),
                    wavelength_nm=float(row["wavelength_nm"]),
                    medium=row["medium"].strip().lower(),
                    source=row["source"],
                    source_url=row["source_url"],
                    source_access_date=row["source_access_date"],
                    source_table_or_query=row["source_table_or_query"],
                    notes=row["notes"],
                )
            )
    return rows


def conventional_rydberg_baseline(path: Path | Traversable | None = None) -> list[dict[str, object]]:
    """Return reference rows with vacuum residuals only where media match."""

    rows: list[dict[str, object]] = []
    for reference in load_reference_lines(path):
        calculated_vacuum_nm = vacuum_transition_wavelength_nm(reference.n_i, reference.n_f)
        comparison_valid = reference.medium == "vacuum"
        residual_nm = calculated_vacuum_nm - reference.wavelength_nm if comparison_valid else None
        relative_error_ppm = (residual_nm / reference.wavelength_nm) * 1_000_000.0 if residual_nm is not None else None
        rows.append(
            {
                "series": reference.series,
                "transition": reference.transition,
                "n_i": reference.n_i,
                "n_f": reference.n_f,
                "raw_reference_wavelength_nm": reference.wavelength_nm,
                "reference_medium": reference.medium,
                "calculated_vacuum_wavelength_nm": calculated_vacuum_nm,
                "comparison_valid": comparison_valid,
                "residual_nm": residual_nm,
                "relative_error_ppm": relative_error_ppm,
                "source": reference.source,
                "source_url": reference.source_url,
                "source_access_date": reference.source_access_date,
            }
        )
    return rows
