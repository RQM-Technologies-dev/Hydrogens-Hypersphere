"""Generate the Hydrogen S^3 Representation diagnostic report."""

from __future__ import annotations

from pathlib import Path

from hydrogen_s3.branching import branching_diagnostics, target_basis
from hydrogen_s3.clebsch_gordan import clebsch_gordan_exact
from hydrogen_s3.reference_spectroscopy import conventional_rydberg_baseline
from hydrogen_s3.spectrum import (
    angular_labels,
    calibrated_spectral_energy_from_K,
    s3_laplacian_eigenvalue,
    shell_dimension,
    shell_number,
    shifted_laplacian_eigenvalue,
)

REPORT_DIR = Path("build/reports")
REPORT_PATH = REPORT_DIR / "hydrogen_s3_report.md"


def _markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def _format_float(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.3e}"
    return str(value)


def build_report() -> str:
    """Return the complete Markdown report."""

    shell_rows = [
        [
            K,
            shell_number(K),
            s3_laplacian_eigenvalue(K),
            shifted_laplacian_eigenvalue(K),
            shell_dimension(K),
            f"{calibrated_spectral_energy_from_K(K):.9f}",
        ]
        for K in range(0, 9)
    ]
    angular_rows = [[K, len(angular_labels(K)), len(target_basis(K))] for K in range(0, 9)]
    diagnostic_rows = []
    for K in range(0, 9):
        diag = branching_diagnostics(K)
        diagnostic_rows.append(
            [
                K,
                diag["dimension"],
                _format_float(diag["unitarity_U_Udagger_error"]),
                _format_float(diag["unitarity_Udagger_U_error"]),
                _format_float(diag["intertwining_Jx_error"]),
                _format_float(diag["intertwining_Jy_error"]),
                _format_float(diag["intertwining_Jz_error"]),
                _format_float(diag["casimir_intertwining_error"]),
            ]
        )
    cg_rows = [
        ["<1/2,1/2;1/2,1/2|1,1>", clebsch_gordan_exact(0.5, 0.5, 0.5, 0.5, 1, 1)],
        ["<1,1;1,0|2,1>", clebsch_gordan_exact(1, 1, 1, 0, 2, 1)],
        ["<1,1;1,-1|0,0>", clebsch_gordan_exact(1, 1, 1, -1, 0, 0)],
    ]
    spectroscopy_rows = []
    for row in conventional_rydberg_baseline():
        spectroscopy_rows.append(
            [
                row["series"],
                row["transition"],
                row["raw_reference_wavelength_nm"],
                row["reference_medium"],
                f"{float(row['calculated_vacuum_wavelength_nm']):.6f}",
                row["comparison_valid"],
                "" if row["residual_nm"] is None else f"{float(row['residual_nm']):.6f}",
            ]
        )

    return (
        "\n\n".join(
            [
                "# Hydrogen S^3 Representation Report",
                (
                    "This generated report covers shell-level representation diagnostics. "
                    "The energy column uses H_spec = -Ry / (-Delta_S3 + 1) as a calibrated "
                    "packaging of the conventional nonrelativistic Rydberg shell law."
                ),
                "## Shell Number and Dimension",
                _markdown_table(["K", "n", "-Delta_S3", "-Delta_S3+1", "dimension", "H_spec eV"], shell_rows),
                "## Angular-Sector Counts",
                _markdown_table(["K", "angular label count", "target basis count"], angular_rows),
                "## Branching Transform Diagnostics",
                _markdown_table(
                    [
                        "K",
                        "dimension",
                        "U Udagger",
                        "Udagger U",
                        "Jx intertwining",
                        "Jy intertwining",
                        "Jz intertwining",
                        "Casimir",
                    ],
                    diagnostic_rows,
                ),
                "## Clebsch-Gordan Convention Checks",
                _markdown_table(["coefficient", "SymPy exact value"], cg_rows),
                "## Appendix: Conventional Rydberg Baseline",
                (
                    "The calculated wavelength is a vacuum value from the supplied Rydberg "
                    "shell law. Residuals are blank when the raw reference medium is not vacuum."
                ),
                _markdown_table(
                    [
                        "series",
                        "transition",
                        "raw reference nm",
                        "medium",
                        "calculated vacuum nm",
                        "comparison valid",
                        "residual nm",
                    ],
                    spectroscopy_rows,
                ),
            ]
        )
        + "\n"
    )


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(build_report(), encoding="utf-8")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
