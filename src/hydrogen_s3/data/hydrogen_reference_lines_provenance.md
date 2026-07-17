# Hydrogen Reference Line Provenance

## Status

The `data/hydrogen_reference_lines.csv` fixture contains a small retained set of
hydrogen reference wavelengths with raw medium labels preserved per row.

## Authoritative Source

- Source: NIST Atomic Spectra Database (ASD), Lines Data for neutral hydrogen (H I)
- URL: https://physics.nist.gov/asd
- Access date for current table values: 2026-05-04

## Query Method Used

1. Open NIST ASD Lines Data from https://physics.nist.gov/asd.
2. Select spectrum `H I`.
3. Use line output with Ritz wavelengths enabled or preferred.
4. Match the transition subset used by this repository using lower and upper principal shells:
   - Lyman: 2->1, 3->1, 4->1, 5->1, 6->1
   - Balmer: 3->2, 4->2, 5->2, 6->2, 7->2
   - Paschen: 4->3, 5->3, 6->3, 7->3, 8->3
5. Record wavelengths in nm together with the wavelength medium exactly as returned by NIST ASD.

## Ritz vs Observed Policy

This dataset uses Ritz wavelengths where NIST provides them.

## Air/Vacuum Handling

- The `medium` column is explicit per row.
- UV Lyman rows are recorded as `vacuum`.
- Optical/IR Balmer and Paschen rows are recorded as `air`.
- The package computes a conventional vacuum wavelength from the supplied Rydberg shell law.
- Residual and relative-error fields are suppressed for non-vacuum reference rows.
- No air-vacuum conversion formula is implemented in this repository.

## Scope Boundary

This fixture supports a conventional Rydberg baseline only. It is not validation
that the S^3 representation derives the Rydberg law or the Coulomb interaction.
