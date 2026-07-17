# Conventional Rydberg Baseline

The repository retains `data/hydrogen_reference_lines.csv` as a small auditable
fixture of hydrogen reference lines. The companion provenance file records the
source, access date, selected transitions, and raw medium labels.

The package computes a conventional vacuum wavelength from the supplied Rydberg
shell law. Because the CSV contains both air and vacuum wavelengths, residuals
are only reported when the reference row is also vacuum. For air rows, the raw
reference wavelength and medium are preserved, but residual and relative-error
fields are left blank.

This avoids comparing a vacuum calculation against an air reference without a
verified conversion. No air-vacuum conversion formula is implemented here.

The baseline is not validation of the S^3 representation. It is a conventional
hydrogen reference calculation using the same supplied shell energy law.
