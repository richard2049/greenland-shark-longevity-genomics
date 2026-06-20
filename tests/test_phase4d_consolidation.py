from pathlib import Path

from greenland_shark_longevity.phase4d_consolidation import (
    build_phase4d_interpretations,
    consolidate_integrated_evidence,
)
from greenland_shark_longevity.schemas import (
    EVIDENCE_COLUMNS,
    PHASE4C_GENE_REVIEW_COLUMNS,
    PHASE4C_LOCUS_REVIEW_COLUMNS,
    PHASE4C_TP53_SUMMARY_COLUMNS,
)
from greenland_shark_longevity.utils import write_tsv


def write_phase4d_fixture_tables(tmp_path: Path) -> tuple[Path, Path, Path]:
    gene_summary = tmp_path / "phase4c_gene_summary.tsv"
    locus_review = tmp_path / "phase4c_locus_review.tsv"
    tp53_summary = tmp_path / "tp53_summary.tsv"

    write_tsv(
        gene_summary,
        [
            {
                "gene_symbol": "FTH1B",
                "mechanism": "Ferroptosis, iron handling, and oxidative stress",
                "reviewed_locus_count": "12",
                "high_coverage_locus_count": "12",
                "domain_supported_locus_count": "12",
                "loci_with_disruption_flags": "1",
                "best_query_coverage": "1",
                "candidate_loci": "scaf1:10-100;scaf2:20-200",
                "phase4c_status": "MULTIPLE_CANDIDATE_LOCI_REQUIRE_MANUAL_REVIEW",
                "artifact_risk": "high",
                "claim_language_guardrail": "test",
                "required_validation": "Manual ferritin-family locus review.",
                "notes": "test",
            },
            {
                "gene_symbol": "H1F0",
                "mechanism": "Chromatin and epigenome-related regulation",
                "reviewed_locus_count": "2",
                "high_coverage_locus_count": "2",
                "domain_supported_locus_count": "2",
                "loci_with_disruption_flags": "0",
                "best_query_coverage": "1",
                "candidate_loci": "scaf3:10-100;scaf4:20-200",
                "phase4c_status": "MULTIPLE_CANDIDATE_LOCI_REQUIRE_MANUAL_REVIEW",
                "artifact_risk": "moderate",
                "claim_language_guardrail": "test",
                "required_validation": "Manual H1-family locus review.",
                "notes": "test",
            },
            {
                "gene_symbol": "RAD51",
                "mechanism": "DNA repair/genome stability",
                "reviewed_locus_count": "2",
                "high_coverage_locus_count": "2",
                "domain_supported_locus_count": "2",
                "loci_with_disruption_flags": "0",
                "best_query_coverage": "1",
                "candidate_loci": "scaf5:10-100;scaf6:20-200",
                "phase4c_status": "MULTIPLE_CANDIDATE_LOCI_REQUIRE_MANUAL_REVIEW",
                "artifact_risk": "moderate",
                "claim_language_guardrail": "test",
                "required_validation": "Manual RAD51-family locus review.",
                "notes": "test",
            },
        ],
        PHASE4C_GENE_REVIEW_COLUMNS,
    )
    write_tsv(
        locus_review,
        [
            {
                "gene_symbol": gene,
                "mechanism": "test",
                "locus_cluster_id": f"{gene}_1",
                "target_seqid": "scaf1",
                "requested_scaffold_aliases": "scaf1",
                "cluster_start": "10",
                "cluster_end": "100",
                "strand_values": "+",
                "supporting_query_count": "1",
                "supporting_query_ids": "q1",
                "supporting_original_protein_ids": "p1",
                "best_query_coverage": "1",
                "best_identity": "0.9",
                "hit_count": "1",
                "high_coverage_hit_count": "1",
                "cds_feature_count_range": "4-4",
                "exon_feature_count_range": "0-0",
                "paf_frameshift_count": "0",
                "paf_stop_count": "0",
                "miniprot_disruption_status": "NO_PARSED_FS_OR_ST_TAG",
                "domain_validation_status": "DOMAIN_SUPPORTED",
                "domain_supported_query_count": "1",
                "best_domain_hits": "PF00000",
                "overlap_cluster_size": "1",
                "overlap_note": "test",
                "candidate_locus_status": "CANDIDATE_LOCUS_SUPPORTED_REQUIRES_REVIEW",
                "artifact_risk": "moderate",
                "claim_language_guardrail": "test",
                "required_validation": "test",
                "notes": "test",
            }
            for gene in ["FTH1B", "H1F0", "RAD51"]
        ],
        PHASE4C_LOCUS_REVIEW_COLUMNS,
    )
    write_tsv(
        tp53_summary,
        [
            {
                "gene_symbol": "TP53",
                "query_id": "TP53__comparator_query__XP_1",
                "query_source_type": "comparator_query_for_unresolved_gene",
                "query_length_aa": "357",
                "search_scope": "phase4_forward_hit_target_scaffolds",
                "chunks_scanned": "1",
                "total_bases_scanned": "260641950",
                "alignment_hit_count": "1",
                "high_coverage_hit_count": "0",
                "partial_hit_count": "1",
                "best_query_coverage": "0.960784",
                "candidate_loci": "scaf247:503768-525127",
                "genome_search_status": "POSSIBLE_DISRUPTED_TP53_ALIGNMENT",
                "candidate_classification": "possible_disrupted_or_divergent_candidate_locus",
                "artifact_risk": "high",
                "claim_language_guardrail": "test",
                "required_validation": "Manual TP53-family review.",
                "notes": "test",
            }
        ],
        PHASE4C_TP53_SUMMARY_COLUMNS,
    )
    return gene_summary, locus_review, tp53_summary


def test_phase4d_keeps_candidate_loci_separate_from_validated_duplication(tmp_path: Path):
    gene_summary, locus_review, tp53_summary = write_phase4d_fixture_tables(tmp_path)

    rows = build_phase4d_interpretations(gene_summary, locus_review, tp53_summary)
    fth1b = next(row for row in rows if row["gene_symbol"] == "FTH1B")

    assert fth1b["evidence_tier"] == "Tier 2"
    assert fth1b["copy_number_interpretation_status"] == "NOT_VALIDATED_DUPLICATION"
    assert "copy-number expansion is not validated" in fth1b["conservative_interpretation"]
    assert fth1b["artifact_risk"] == "high"


def test_phase4d_tp53_remains_artifact_uncertain_without_gene_state_claim(tmp_path: Path):
    gene_summary, locus_review, tp53_summary = write_phase4d_fixture_tables(tmp_path)

    rows = build_phase4d_interpretations(gene_summary, locus_review, tp53_summary)
    tp53 = next(row for row in rows if row["gene_symbol"] == "TP53")

    assert tp53["evidence_tier"] == "Artifact/uncertain"
    assert tp53["phase4d_status"] == "TP53_CANDIDATE_ALIGNMENT_UNCERTAIN_REQUIRES_REVIEW"
    lowered = tp53["conservative_interpretation"].lower()
    assert "inactivated" not in lowered
    assert "absent" not in lowered
    assert "establish gene state" in lowered


def test_phase4d_replaces_earlier_high_priority_evidence_rows(tmp_path: Path):
    gene_summary, locus_review, tp53_summary = write_phase4d_fixture_tables(tmp_path)
    base_evidence = tmp_path / "base_evidence.tsv"
    write_tsv(
        base_evidence,
        [
            {
                "mechanism": "DNA repair/genome stability",
                "gene_or_pathway": "RAD51",
                "evidence_type": "candidate_presence_or_not_assessed",
                "evidence_tier": "Artifact/uncertain",
                "resources_supporting": "NOT_ASSESSED",
                "artifact_risk": "not_assessable",
                "biological_interpretation": "Do not claim biologically; evidence is incomplete or artifact-prone.",
                "relevance_to_aging_longevity": "NOT_ASSESSED",
                "translational_relevance": "NOT_ASSESSED",
                "required_validation": "REQUIRES_VALIDATION",
                "claim_language_guardrail": "test",
            },
            {
                "mechanism": "DNA repair/genome stability",
                "gene_or_pathway": "BRCA1",
                "evidence_type": "candidate_presence_or_not_assessed",
                "evidence_tier": "Tier 2",
                "resources_supporting": "SMIC_TOKYO_GENOME_2025",
                "artifact_risk": "moderate",
                "biological_interpretation": "Plausible but incomplete computational evidence.",
                "relevance_to_aging_longevity": "Hypothesis-generating; requires validation before biological interpretation.",
                "translational_relevance": "NOT_ASSESSED",
                "required_validation": "REQUIRES_VALIDATION",
                "claim_language_guardrail": "test",
            },
        ],
        EVIDENCE_COLUMNS,
    )

    interpretations = build_phase4d_interpretations(gene_summary, locus_review, tp53_summary)
    consolidated = consolidate_integrated_evidence(base_evidence, interpretations)

    rad51_rows = [row for row in consolidated if row["gene_or_pathway"] == "RAD51"]
    assert len(rad51_rows) == 1
    assert rad51_rows[0]["evidence_type"] == "phase4d_candidate_locus_review"
    assert any(row["gene_or_pathway"] == "BRCA1" for row in consolidated)
