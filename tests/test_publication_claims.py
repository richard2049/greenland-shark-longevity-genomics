from pathlib import Path

from greenland_shark_longevity.publication_claims import audit_publication_claims


def test_publication_claim_audit_blocks_unreproduced_pnas_mechanisms():
    rows = audit_publication_claims(
        Path("config/publication_claims.yaml"),
        Path("data/metadata/reference_file_inventory.tsv"),
    )
    by_id = {row["claim_id"]: row for row in rows}
    assert by_id["PNAS_YANG_2026_RESOURCE_ASSEMBLY"]["evidence_tier"] == "Resource-quality observation"
    assert by_id["PNAS_YANG_2026_LINKER_HISTONE"]["evidence_tier"] == "NOT_ASSESSED"
    assert "orthology" in by_id["PNAS_YANG_2026_LINKER_HISTONE"]["blockers"]
    assert "genome_wide_protein_set" not in by_id["PNAS_YANG_2026_LINKER_HISTONE"]["blockers"]
    assert "separable_loci" in by_id["PNAS_YANG_2026_FTH1B_COPY_NUMBER"]["blockers"]
    assert "annotation_coordinates" not in by_id["PNAS_YANG_2026_FTH1B_COPY_NUMBER"]["blockers"]
