from greenland_shark_longevity.phase8b_final_evidence_scoring import integrate_phase8b


def evidence_row(gene: str, tier: str = "Tier 2", artifact: str = "low", evidence_type: str = "candidate_presence_or_not_assessed") -> dict[str, str]:
    return {
        "mechanism": "DNA repair/genome stability",
        "gene_or_pathway": gene,
        "evidence_type": evidence_type,
        "evidence_tier": tier,
        "resources_supporting": "SMIC_TOKYO_GENOME_2025",
        "artifact_risk": artifact,
        "biological_interpretation": "Plausible but incomplete computational evidence.",
        "relevance_to_aging_longevity": "Hypothesis-generating; requires validation before biological interpretation.",
        "translational_relevance": "NOT_ASSESSED",
        "required_validation": "REQUIRES_VALIDATION",
        "claim_language_guardrail": "Avoid activated/inactivated/absent/causal wording unless direct evidence supports it.",
    }


def expression_audit_row(gene: str, support_used: str = "CAUTIOUS_TISSUE_SPECIFIC_SUPPORT") -> dict[str, str]:
    return {
        "gene_symbol": gene,
        "base_evidence_tier": "Tier 2",
        "expression_evidence_tier": "Tier 2",
        "integrated_evidence_tier": "Tier 2",
        "tier_change": "UNCHANGED",
        "phase7e_expression_support_level": "CAUTIOUS_RETINA_EXPRESSION_SUPPORT_FOR_PHASE8",
        "expression_interpretation_status": "RETINA_DETECTED_FIRST_PASS_WITH_QC_CAVEATS",
        "integration_action": "APPEND_CAUTION_RETINA_EXPRESSION_SUPPORT_NO_TIER_CHANGE",
        "expression_support_used": support_used,
        "artifact_risk_update": "low_mapping_rate_caveat",
        "required_validation_update": "Review expression reference.",
        "supporting_files": "results/rnaseq/phase7e_candidate_expression_hardened.tsv",
    }


def phase4e_row(gene: str, copy_status: str = "NOT_VALIDATED_DUPLICATION", artifact: str = "moderate") -> dict[str, str]:
    return {
        "gene_symbol": gene,
        "mechanism": "DNA repair/genome stability",
        "reviewed_locus_count": "2",
        "focal_annotation_overlap_count": "2",
        "product_consistent_annotation_overlap_count": "2",
        "domain_supported_locus_count": "2",
        "no_disruption_locus_count": "2",
        "disruption_flag_locus_count": "0",
        "distinct_annotation_seqids": "scaffold_1",
        "min_inter_locus_distance_bp": "1000",
        "clustered_same_scaffold": "no",
        "phase4e_hardened_status": "HARDENED_CANDIDATE_LOCI_REQUIRE_CROSS_RESOURCE_VALIDATION",
        "evidence_tier_recommendation": "Tier 2",
        "artifact_risk": artifact,
        "copy_number_interpretation_status": copy_status,
        "conservative_interpretation": "Candidate loci require review.",
        "next_action": "Review loci.",
        "required_validation": "Manual locus review and cross-resource validation.",
        "claim_language_guardrail": "Do not claim validated duplication.",
        "supporting_files": "results/rescue/phase4e_gene_hardened_summary.tsv",
    }


def phase5c_row(gene: str, artifact: str = "high_repeat_context_artifact_risk") -> dict[str, str]:
    return {
        "gene_symbol": gene,
        "mechanism": "DNA repair/genome stability",
        "locus_count": "2",
        "annotation_seqids": "scaffold_1",
        "gff_loci_with_direct_overlap": "2",
        "repeatmasker_out_loci_with_direct_overlap": "2",
        "max_repeatmasker_out_overlap_fraction": "0.8",
        "max_window_repeat_density": "0.8",
        "repeatmasker_out_classes": "LINE",
        "repeatmasker_out_families": "CR1",
        "gff_out_concordance_summary": "BOTH_GFF_AND_OUT_DIRECT_OVERLAP",
        "phase5c_qc_status": "REPEAT_CONTEXT_QC_SUPPORTED_BY_GFF_AND_OUT_ARTIFACT_ONLY",
        "artifact_risk_modifier": artifact,
        "conservative_interpretation": "Repeat context is artifact evidence only.",
        "required_validation": "Review repeat intervals.",
        "claim_language_guardrail": "Do not interpret repeats as mechanism.",
        "supporting_files": "results/repeats/phase5c_gene_repeat_qc_summary.tsv",
    }


def phase6_row(gene: str) -> dict[str, str]:
    return {
        "mechanism": "Telomere-related biology",
        "gene_symbol": gene,
        "synonyms": "NOT_ASSESSED",
        "panel_status": "CURATED_PHASE6_TELOMERE_PANEL",
        "copy_number_mapping_status": "ANNOTATION_SYMBOL_MATCH_ORTHOGROUP_MAPPED",
        "copy_count": "1",
        "orthogroup_id": "OG1",
        "resources_supporting": "SMIC_TOKYO_GENOME_2025",
        "duplication_artifact_risk": "low",
        "integrated_evidence_tier": "Tier 2",
        "integrated_artifact_risk": "low",
        "readiness_status": "ORTHOLOGY_EVIDENCE_READY_FOR_DOMAIN_AND_CONTEXT_VALIDATION",
        "conservative_interpretation": "Readiness only.",
        "required_validation": "Validate telomere gene context.",
        "claim_language_guardrail": "Do not infer telomere length.",
        "supporting_files": "results/telomere/phase6_telomere_gene_audit.tsv",
    }


def test_phase8b_downgrades_unsupported_tier1_to_tier2():
    final_rows, audit_rows, _ = integrate_phase8b(
        [evidence_row("GENE1", tier="Tier 1", artifact="high")],
        [expression_audit_row("GENE1")],
        [],
        [],
        [],
    )

    assert final_rows[0]["evidence_tier"] == "Tier 2"
    assert audit_rows[0]["tier_change"] == "Tier 1->Tier 2"
    assert audit_rows[0]["phase8b_audit_status"] == "TIER_DOWNGRADED_BY_PHASE8B_RULES"


def test_phase8b_retains_high_risk_tier2_as_artifact_prone_plausible_lead():
    final_rows, audit_rows, _ = integrate_phase8b(
        [evidence_row("RAD51", tier="Tier 2", artifact="moderate;high_repeat_context_artifact_risk")],
        [expression_audit_row("RAD51", support_used="NO_POSITIVE_EXPRESSION_SUPPORT_USED")],
        [phase4e_row("RAD51")],
        [phase5c_row("RAD51")],
        [],
    )

    assert final_rows[0]["evidence_tier"] == "Tier 2"
    assert audit_rows[0]["phase8b_reporting_class"] == "PLAUSIBLE_LEAD_ARTIFACT_PRONE"
    assert "high_repeat_context_artifact_risk" in audit_rows[0]["major_caveats"]
    assert "phase8b_final_tier_audit" in final_rows[0]["evidence_type"]


def test_phase8b_retains_artifact_uncertain_as_do_not_claim():
    final_rows, audit_rows, _ = integrate_phase8b(
        [evidence_row("TP53", tier="Artifact/uncertain", artifact="high")],
        [expression_audit_row("TP53", support_used="NO_POSITIVE_EXPRESSION_SUPPORT_USED")],
        [],
        [],
        [],
    )

    assert final_rows[0]["evidence_tier"] == "Artifact/uncertain"
    assert audit_rows[0]["phase8b_reporting_class"] == "ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY"


def test_phase8b_mechanism_summary_counts_expression_and_telomere_context():
    rows = [evidence_row("TERT", tier="Tier 2", artifact="low")]
    rows[0]["mechanism"] = "Telomere-related biology"
    final_rows, audit_rows, summary_rows = integrate_phase8b(
        rows,
        [expression_audit_row("TERT")],
        [],
        [],
        [phase6_row("TERT")],
    )

    assert final_rows[0]["evidence_tier"] == "Tier 2"
    assert audit_rows[0]["telomere_component"] == "ORTHOLOGY_EVIDENCE_READY_FOR_DOMAIN_AND_CONTEXT_VALIDATION"
    assert summary_rows[0]["expression_supported_count"] == "1"
