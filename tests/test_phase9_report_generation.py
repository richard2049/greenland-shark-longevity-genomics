from pathlib import Path

from greenland_shark_longevity.phase9_report_generation import generate_phase9_outputs
from greenland_shark_longevity.schemas import PHASE9_FIGURE_MANIFEST_COLUMNS, PHASE9_KEY_FINDINGS_COLUMNS
from greenland_shark_longevity.utils import read_tsv


def audit_row(
    gene_symbol: str,
    mechanism: str,
    tier: str,
    reporting_class: str,
    artifact_risk_level: str,
    major_caveats: str = "standard_validation_required",
) -> dict[str, str]:
    return {
        "gene_symbol": gene_symbol,
        "mechanism": mechanism,
        "current_evidence_tier": tier,
        "phase8b_final_evidence_tier": tier,
        "tier_change": "UNCHANGED",
        "phase8b_audit_status": "TIER_RETAINED_PENDING_VALIDATION",
        "phase8b_reporting_class": reporting_class,
        "evidence_components": "ORTHOLOGY_OR_CANDIDATE_MAPPING_PRESENT",
        "orthology_component": "ORTHOLOGY_OR_CANDIDATE_MAPPING_PRESENT",
        "locus_component": "LOCUS_NOT_TARGETED_BY_PHASE4E",
        "repeat_context_component": "REPEAT_CONTEXT_NOT_TARGETED_BY_PHASE5C",
        "telomere_component": "NOT_TELOMERE_PANEL_CANDIDATE",
        "expression_component": "CAUTIOUS_TISSUE_SPECIFIC_SUPPORT",
        "artifact_risk_level": artifact_risk_level,
        "major_caveats": major_caveats,
        "required_validation": "REQUIRES_VALIDATION",
        "supporting_files": "results/evidence/phase8b_tier_audit.tsv",
        "claim_language_guardrail": "Use cautious candidate language only.",
    }


def mechanism_row(
    mechanism: str,
    candidate_count: int,
    tier2_count: int,
    artifact_uncertain_count: int,
    artifact_prone_count: int,
    expression_supported_count: int,
    repeat_context_caveat_count: int,
    locus_review_required_count: int,
) -> dict[str, str]:
    return {
        "mechanism": mechanism,
        "candidate_count": str(candidate_count),
        "tier1_count": "0",
        "tier2_count": str(tier2_count),
        "tier3_count": "0",
        "artifact_uncertain_count": str(artifact_uncertain_count),
        "robust_ready_count": "0",
        "plausible_lead_count": str(tier2_count),
        "exploratory_count": "0",
        "artifact_prone_count": str(artifact_prone_count),
        "expression_supported_count": str(expression_supported_count),
        "repeat_context_caveat_count": str(repeat_context_caveat_count),
        "locus_review_required_count": str(locus_review_required_count),
        "required_validation": "Use row-level audit entries.",
    }


def phase9_paths(tmp_path: Path) -> dict[str, Path]:
    return {
        "report_output": tmp_path / "reports" / "final" / "report.md",
        "figure_manifest_output": tmp_path / "results" / "reporting" / "phase9_figure_manifest.tsv",
        "key_findings_output": tmp_path / "results" / "reporting" / "phase9_key_findings.tsv",
        "tier_data_output": tmp_path / "reports" / "figures" / "data" / "phase9_evidence_tier_summary.tsv",
        "class_data_output": tmp_path / "reports" / "figures" / "data" / "phase9_reporting_class_summary.tsv",
        "mechanism_matrix_data_output": tmp_path / "reports" / "figures" / "data" / "phase9_mechanism_evidence_matrix.tsv",
        "artifact_context_data_output": tmp_path / "reports" / "figures" / "data" / "phase9_artifact_context_summary.tsv",
        "tier_figure_output": tmp_path / "reports" / "figures" / "phase9_evidence_tier_summary.svg",
        "class_figure_output": tmp_path / "reports" / "figures" / "phase9_reporting_class_summary.svg",
        "mechanism_matrix_figure_output": tmp_path / "reports" / "figures" / "phase9_mechanism_evidence_matrix.svg",
        "artifact_context_figure_output": tmp_path / "reports" / "figures" / "phase9_artifact_context_summary.svg",
    }


def test_phase9_generates_manifest_figures_and_report(tmp_path: Path) -> None:
    audit_rows = [
        audit_row("BRCA1", "DNA repair/genome stability", "Tier 2", "PLAUSIBLE_LEAD_REQUIRES_VALIDATION", "low"),
        audit_row("FTH1B", "Ferroptosis, iron handling, and oxidative stress", "Tier 2", "PLAUSIBLE_LEAD_ARTIFACT_PRONE", "high", "high_repeat_context_artifact_risk"),
        audit_row("TP53", "p53 pathway", "Artifact/uncertain", "ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY", "high"),
    ]
    summary_rows = [
        mechanism_row("DNA repair/genome stability", 1, 1, 0, 0, 1, 0, 0),
        mechanism_row("Ferroptosis, iron handling, and oxidative stress", 1, 1, 0, 1, 0, 1, 1),
        mechanism_row("p53 pathway", 1, 0, 1, 1, 0, 1, 1),
    ]
    paths = phase9_paths(tmp_path)

    figure_manifest, key_findings = generate_phase9_outputs(audit_rows, summary_rows, **paths)

    assert len(figure_manifest) == 4
    assert len(key_findings) == 5
    for path in paths.values():
        assert path.exists(), path
    manifest_rows = read_tsv(paths["figure_manifest_output"], PHASE9_FIGURE_MANIFEST_COLUMNS)
    assert {row["status"] for row in manifest_rows} == {"GENERATED"}
    key_rows = read_tsv(paths["key_findings_output"], PHASE9_KEY_FINDINGS_COLUMNS)
    assert "No current candidate reaches Tier 1" in key_rows[0]["answer"]
    assert "FTH1B" in key_rows[2]["answer"]
    assert "TP53" in key_rows[3]["answer"]
    mechanism_svg = paths["mechanism_matrix_figure_output"].read_text(encoding="utf-8")
    assert "<svg" in mechanism_svg
    assert "<title" in mechanism_svg
    assert "Candidate count" in mechanism_svg
    assert "Data source:" in mechanism_svg
    assert "Guardrail:" in mechanism_svg
    assert "Context counters" in mechanism_svg


def test_phase9_report_keeps_interpretation_conservative(tmp_path: Path) -> None:
    audit_rows = [
        audit_row("RAD51", "DNA repair/genome stability", "Tier 2", "PLAUSIBLE_LEAD_ARTIFACT_PRONE", "high", "copy_number_not_validated"),
        audit_row("TP53", "p53 pathway", "Artifact/uncertain", "ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY", "high"),
    ]
    summary_rows = [
        mechanism_row("DNA repair/genome stability", 1, 1, 0, 1, 0, 1, 1),
        mechanism_row("p53 pathway", 1, 0, 1, 1, 0, 1, 1),
    ]
    paths = phase9_paths(tmp_path)

    generate_phase9_outputs(audit_rows, summary_rows, **paths)

    report = paths["report_output"].read_text(encoding="utf-8")
    assert "does not change evidence tiers" in report
    assert "Repeat context is artifact/context evidence only" in report
    assert "TP53` remains artifact/uncertain" in report
    assert "functional advantage" in report
