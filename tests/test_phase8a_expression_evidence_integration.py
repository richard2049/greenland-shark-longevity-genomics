from greenland_shark_longevity.phase8a_expression_evidence_integration import (
    CAUTIOUS_SUPPORT,
    DEFER_SUPPORT,
    LIMITED_SUPPORT,
    build_expression_evidence_row,
    integrate_phase8a_expression_support,
)


def base_evidence_row(gene: str, tier: str = "Tier 2") -> dict[str, str]:
    return {
        "mechanism": "DNA repair/genome stability",
        "gene_or_pathway": gene,
        "evidence_type": "candidate_presence_or_not_assessed",
        "evidence_tier": tier,
        "resources_supporting": "SMIC_TOKYO_GENOME_2025",
        "artifact_risk": "low",
        "biological_interpretation": "Plausible but incomplete computational evidence.",
        "relevance_to_aging_longevity": "Hypothesis-generating; requires validation before biological interpretation.",
        "translational_relevance": "NOT_ASSESSED",
        "required_validation": "REQUIRES_VALIDATION",
        "claim_language_guardrail": "Avoid activated/inactivated/absent/causal wording unless direct evidence supports it.",
    }


def expression_row(gene: str, support_level: str, status: str | None = None) -> dict[str, str]:
    return {
        "mechanism": "DNA repair/genome stability",
        "gene_symbol": gene,
        "tissue": "retina",
        "quantified_run_count": "3",
        "detected_run_count": "3" if support_level == CAUTIOUS_SUPPORT else "1",
        "not_detected_run_count": "0",
        "locus_review_run_count": "0" if support_level != DEFER_SUPPORT else "3",
        "detected_runs": "SRR1;SRR2;SRR3",
        "median_tpm": "5.0",
        "min_tpm": "1.0",
        "max_tpm": "9.0",
        "median_numreads": "100.0",
        "phase7c_quantification_readiness": "REFERENCE_READY_FOR_FIRST_PASS_QUANTIFICATION",
        "reference_mapping_status": "TRANSCRIPT_REFERENCE_PRESENT",
        "reference_ambiguity_status": "low_reference_mapping_risk",
        "matched_transcript_count": "1",
        "orthology_mapping_status": "ANNOTATION_SYMBOL_MATCH_ORTHOGROUP_MAPPED",
        "duplication_artifact_risk": "low",
        "phase4e_hardened_status": "NOT_ASSESSED",
        "phase4e_artifact_risk": "NOT_ASSESSED",
        "phase5c_artifact_risk_modifier": "NOT_ASSESSED",
        "low_mapping_run_count": "2",
        "very_low_mapping_run_count": "1",
        "min_run_mapping_rate": "14.0",
        "expression_interpretation_status": status or "RETINA_DETECTED_FIRST_PASS_WITH_QC_CAVEATS",
        "phase7e_expression_support_level": support_level,
        "artifact_risk": "low_reference_mapping_risk;low_mapping_rate_caveat",
        "conservative_interpretation": "Fixture interpretation.",
        "required_validation": "Review expression reference.",
        "claim_language_guardrail": "Use Phase 7e only for retina-specific expression support.",
        "supporting_files": "results/rnaseq/phase7e_candidate_expression_hardened.tsv",
    }


def test_cautious_expression_support_does_not_upgrade_integrated_tier():
    base_rows = [base_evidence_row("ERCC4", "Tier 2")]
    expression_rows = [expression_row("ERCC4", CAUTIOUS_SUPPORT)]

    expression_evidence, audit_rows, integrated = integrate_phase8a_expression_support(base_rows, expression_rows)

    assert expression_evidence[0]["evidence_tier"] == "Tier 2"
    assert integrated[0]["evidence_tier"] == "Tier 2"
    assert audit_rows[0]["tier_change"] == "UNCHANGED"
    assert "phase8a_retina_expression_support_cautious" in integrated[0]["evidence_type"]
    assert "retina-specific expression detection" in integrated[0]["biological_interpretation"]


def test_expression_support_cannot_rescue_artifact_uncertain_base_tier():
    base_rows = [base_evidence_row("TP53", "Artifact/uncertain")]
    expression_rows = [expression_row("TP53", CAUTIOUS_SUPPORT)]

    _, audit_rows, integrated = integrate_phase8a_expression_support(base_rows, expression_rows)

    assert integrated[0]["evidence_tier"] == "Artifact/uncertain"
    assert audit_rows[0]["base_evidence_tier"] == "Artifact/uncertain"
    assert audit_rows[0]["integrated_evidence_tier"] == "Artifact/uncertain"


def test_deferred_locus_ambiguous_expression_is_recorded_as_caveat_not_positive_support():
    base_rows = [base_evidence_row("RAD51", "Tier 2")]
    expression_rows = [
        expression_row(
            "RAD51",
            DEFER_SUPPORT,
            status="LOCUS_REVIEW_REQUIRED_BEFORE_EXPRESSION_SUPPORT",
        )
    ]

    expression_evidence, audit_rows, integrated = integrate_phase8a_expression_support(base_rows, expression_rows)

    assert expression_evidence[0]["evidence_tier"] == "Artifact/uncertain"
    assert audit_rows[0]["expression_support_used"] == "NO_POSITIVE_EXPRESSION_SUPPORT_USED"
    assert integrated[0]["evidence_tier"] == "Tier 2"
    assert "expression_support_deferred_locus_or_reference_ambiguity" in integrated[0]["artifact_risk"]


def test_limited_expression_signal_is_tier3_expression_evidence_only():
    row = expression_row(
        "GENE1",
        LIMITED_SUPPORT,
        status="RETINA_DETECTED_IN_LIMITED_RUNS_EXPLORATORY",
    )

    evidence = build_expression_evidence_row(row)

    assert evidence["evidence_tier"] == "Tier 3"
    assert evidence["evidence_type"] == "phase8a_retina_expression_support_limited_exploratory"
