from pathlib import Path

from greenland_shark_longevity.phase5_evidence_integration import integrate_phase5_repeat_context
from greenland_shark_longevity.schemas import EVIDENCE_COLUMNS, PHASE5_GENE_REPEAT_CONTEXT_COLUMNS
from greenland_shark_longevity.utils import write_tsv


def evidence_row(gene: str, tier: str, artifact_risk: str = "moderate") -> dict[str, str]:
    return {
        "mechanism": "DNA repair/genome stability" if gene == "RAD51" else "p53 pathway",
        "gene_or_pathway": gene,
        "evidence_type": "phase4e_manual_locus_hardening",
        "evidence_tier": tier,
        "resources_supporting": "SMIC_TOKYO_GENOME_2025;Pfam-A_38.1",
        "artifact_risk": artifact_risk,
        "biological_interpretation": f"{gene} is a candidate locus requiring validation.",
        "relevance_to_aging_longevity": "Hypothesis-generating only; requires cross-resource and functional validation."
        if tier == "Tier 2"
        else "NOT_ASSESSED",
        "translational_relevance": "NOT_ASSESSED",
        "required_validation": "Manual locus review before biological interpretation.",
        "claim_language_guardrail": "Do not claim validated duplication, function, adaptation, activation, inactivation, loss, causation, or translational relevance.",
    }


def repeat_summary_row(gene: str, overlap_loci: str, max_fraction: str) -> dict[str, str]:
    return {
        "gene_symbol": gene,
        "mechanism": "DNA repair/genome stability" if gene == "RAD51" else "p53 pathway",
        "locus_count": "2" if gene == "RAD51" else "1",
        "annotation_seqids": "scaffold_12" if gene == "RAD51" else "scaffold_247",
        "cluster_start": "100",
        "cluster_end": "500",
        "repeat_annotation_status": "REPEAT_ANNOTATION_AVAILABLE",
        "loci_with_repeat_overlap": overlap_loci,
        "total_repeat_overlap_bp": "1000",
        "max_locus_repeat_overlap_fraction": max_fraction,
        "repeat_classes": "dispersed_repeat",
        "repeat_families": "NONE_IN_WINDOW",
        "phase5_context_status": "LOCAL_REPEAT_CONTEXT_RECORDED_ARTIFACT_CONTEXT_ONLY",
        "artifact_risk": "context_recorded_not_interpreted",
        "conservative_interpretation": "Local repeat annotations were intersected with candidate loci as artifact/context evidence only.",
        "required_validation": "Confirm repeat annotation provenance and candidate-locus coordinates.",
        "claim_language_guardrail": "Do not interpret Phase 5 local repeat context as biological mechanism, validated duplication, adaptation, pathway activity, or longevity evidence.",
        "supporting_files": "results/repeats/phase5_candidate_locus_repeat_context.tsv;results/repeats/phase5_repeat_features.tsv",
    }


def test_phase5_repeat_context_never_upgrades_evidence_tier(tmp_path: Path):
    base = tmp_path / "base.tsv"
    phase4e = tmp_path / "phase4e.tsv"
    summary = tmp_path / "summary.tsv"
    write_tsv(base, [evidence_row("TP53", "Artifact/uncertain", "high")], EVIDENCE_COLUMNS)
    write_tsv(phase4e, [evidence_row("TP53", "Artifact/uncertain", "high")], EVIDENCE_COLUMNS)
    write_tsv(summary, [repeat_summary_row("TP53", "1", "0.766901")], PHASE5_GENE_REPEAT_CONTEXT_COLUMNS)

    phase5_rows, integrated = integrate_phase5_repeat_context(base, phase4e, summary)

    assert phase5_rows[0]["gene_or_pathway"] == "TP53"
    assert phase5_rows[0]["evidence_tier"] == "Artifact/uncertain"
    assert "high_repeat_context_artifact_risk" in phase5_rows[0]["artifact_risk"]
    assert integrated[0]["gene_or_pathway"] == "TP53"
    assert integrated[0]["evidence_tier"] == "Artifact/uncertain"


def test_phase5_repeat_context_replaces_existing_gene_row_and_retains_other_genes(tmp_path: Path):
    base = tmp_path / "base.tsv"
    phase4e = tmp_path / "phase4e.tsv"
    summary = tmp_path / "summary.tsv"
    brca1 = evidence_row("BRCA1", "Tier 2", "low")
    brca1["mechanism"] = "DNA repair/genome stability"
    write_tsv(base, [evidence_row("RAD51", "Tier 2"), brca1], EVIDENCE_COLUMNS)
    write_tsv(phase4e, [evidence_row("RAD51", "Tier 2")], EVIDENCE_COLUMNS)
    write_tsv(summary, [repeat_summary_row("RAD51", "2", "0.829597")], PHASE5_GENE_REPEAT_CONTEXT_COLUMNS)

    phase5_rows, integrated = integrate_phase5_repeat_context(base, phase4e, summary)

    rad51_rows = [row for row in integrated if row["gene_or_pathway"] == "RAD51"]
    assert len(rad51_rows) == 1
    assert "phase5_repeat_context_artifact_risk" in rad51_rows[0]["evidence_type"]
    assert "high_repeat_context_artifact_risk" in rad51_rows[0]["artifact_risk"]
    assert rad51_rows[0]["evidence_tier"] == "Tier 2"
    assert any(row["gene_or_pathway"] == "BRCA1" for row in integrated)
    assert phase5_rows[0]["resources_supporting"].endswith("RepeatModeler_2.0.8;RepeatMasker_4.2.3")
