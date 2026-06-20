from pathlib import Path

from greenland_shark_longevity.evidence import score_audit_row, score_evidence
from greenland_shark_longevity.orthology import integrate_candidates


def test_demo_candidate_integration_outputs_rows():
    copy_rows, audit_rows = integrate_candidates(
        Path("config/candidate_panels.yaml"),
        Path("data/demo/orthology/demo_orthogroups.tsv"),
        Path("data/demo/annotations/demo_gene_coordinates.tsv"),
        Path("data/demo/annotations/demo_domains.tsv"),
        Path("data/demo/transcriptomics/demo_expression_support.tsv"),
        demo_only=True,
    )
    tp53 = next(row for row in copy_rows if row["gene_symbol"] == "TP53")
    assert tp53["copy_count"] == 2
    assert tp53["demo_only"] == "True"
    ercc4 = next(row for row in audit_rows if row["gene_symbol"] == "ERCC4")
    assert ercc4["fragmentation_risk"] == "yes"


def test_demo_evidence_is_not_tier1(tmp_path):
    _copy_rows, audit_rows = integrate_candidates(
        Path("config/candidate_panels.yaml"),
        Path("data/demo/orthology/demo_orthogroups.tsv"),
        Path("data/demo/annotations/demo_gene_coordinates.tsv"),
        Path("data/demo/annotations/demo_domains.tsv"),
        Path("data/demo/transcriptomics/demo_expression_support.tsv"),
        demo_only=True,
    )
    for row in audit_rows:
        assert score_audit_row(row)["evidence_tier"] != "Tier 1"


def test_tier1_requires_orthology_domain_and_context():
    base = {
        "mechanism": "DNA repair/genome stability",
        "gene_symbol": "RAD51",
        "copy_count": "2",
        "orthology_support": "yes",
        "domain_integrity": "complete",
        "separable_loci": "yes",
        "coordinate_support": "yes",
        "cross_resource_support": "no",
        "artifact_risk": "low",
        "resources_supporting": "resource_a",
        "required_validation": "REQUIRES_VALIDATION",
        "demo_only": "false",
    }
    assert score_audit_row(base)["evidence_tier"] == "Tier 1"
    missing_domain = dict(base, domain_integrity="NOT_ASSESSED")
    assert score_audit_row(missing_domain)["evidence_tier"] != "Tier 1"

