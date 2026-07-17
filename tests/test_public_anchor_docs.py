from pathlib import Path


def test_anchor_note_exists_and_contains_required_phrases():
    path = Path("notes/s3_scalar_harmonic_shell_architecture.md")
    assert path.exists()
    text = path.read_text(encoding="utf-8")
    assert "Hydrogen bound-state shell architecture" in text
    assert "scalar harmonics on S" in text
    assert ("Nhat" in text) or ("\\hat N" in text)
    assert "Honesty boundary" in text


def test_claims_matrix_contains_canon_v2_hydrogen_boundary():
    text = Path("docs/claims_matrix.md").read_text(encoding="utf-8")
    assert "Hydrogen shell labels organized on S³" in text
    assert "Standard-compatible" in text
    assert "not Fock's derivation" in text


def test_readme_contains_central_public_claim_section():
    text = Path("README.md").read_text(encoding="utf-8")
    assert "Central public claim" in text
