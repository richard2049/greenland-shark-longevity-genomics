"""Phase 8b final integrated evidence scoring and tier audit."""

from __future__ import annotations

import argparse
import logging
from collections import Counter, defaultdict
from pathlib import Path

from .evidence import validate_guardrail_language
from .phase8a_expression_evidence_integration import append_unique_semicolon, join_sentences
from .schemas import (
    EVIDENCE_COLUMNS,
    PHASE4E_GENE_HARDENING_COLUMNS,
    PHASE5C_GENE_QC_COLUMNS,
    PHASE6_TELOMERE_GENE_AUDIT_COLUMNS,
    PHASE8A_INTEGRATION_AUDIT_COLUMNS,
    PHASE8B_MECHANISM_SUMMARY_COLUMNS,
    PHASE8B_TIER_AUDIT_COLUMNS,
)
from .utils import NOT_ASSESSED, ensure_parent, read_tsv, read_yaml, write_tsv

LOGGER = logging.getLogger(__name__)

TIER_ORDER = {"Artifact/uncertain": 0, "Tier 3": 1, "Tier 2": 2, "Tier 1": 3}
CAUTIOUS_EXPRESSION = "CAUTIOUS_TISSUE_SPECIFIC_SUPPORT"


def index_by(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in rows}


def value_has(value: str | None, *needles: str) -> bool:
    text = (value or "").lower()
    return any(needle.lower() in text for needle in needles)


def artifact_risk_level(row: dict[str, str], phase5c_row: dict[str, str] | None, expression_row: dict[str, str] | None) -> str:
    fields = [
        row.get("artifact_risk", ""),
        phase5c_row.get("artifact_risk_modifier", "") if phase5c_row else "",
        expression_row.get("artifact_risk_update", "") if expression_row else "",
    ]
    text = ";".join(fields).lower()
    if row["evidence_tier"] == "Artifact/uncertain" or "high" in text or "not_assessable" in text:
        return "high"
    if "moderate" in text or "deferred" in text or "not_interpretable" in text:
        return "moderate"
    if "low" in text:
        return "low"
    return "not_assessed"


def orthology_component(row: dict[str, str]) -> str:
    evidence_type = row.get("evidence_type", "")
    if "candidate_presence_or_not_assessed" in evidence_type or "candidate_copy_number" in evidence_type:
        return "ORTHOLOGY_OR_CANDIDATE_MAPPING_PRESENT"
    if row["evidence_tier"] == "Artifact/uncertain":
        return "ORTHOLOGY_UNCERTAIN_OR_NOT_ASSESSABLE"
    return "ORTHOLOGY_NOT_EXPLICITLY_RECORDED"


def locus_component(phase4e_row: dict[str, str] | None) -> str:
    if phase4e_row is None:
        return "LOCUS_NOT_TARGETED_BY_PHASE4E"
    return phase4e_row.get("phase4e_hardened_status", "LOCUS_REVIEW_STATUS_NOT_ASSESSED")


def repeat_component(phase5c_row: dict[str, str] | None) -> str:
    if phase5c_row is None:
        return "REPEAT_CONTEXT_NOT_TARGETED_BY_PHASE5C"
    return phase5c_row.get("phase5c_qc_status", "REPEAT_CONTEXT_STATUS_NOT_ASSESSED")


def telomere_component(phase6_row: dict[str, str] | None) -> str:
    if phase6_row is None:
        return "NOT_TELOMERE_PANEL_CANDIDATE"
    return phase6_row.get("readiness_status", "TELOMERE_READINESS_NOT_ASSESSED")


def expression_component(expression_row: dict[str, str] | None) -> str:
    if expression_row is None:
        return "EXPRESSION_NOT_REVIEWED_BY_PHASE8A"
    return expression_row.get("expression_support_used", "EXPRESSION_SUPPORT_NOT_ASSESSED")


def evidence_components(
    row: dict[str, str],
    phase4e_row: dict[str, str] | None,
    phase5c_row: dict[str, str] | None,
    phase6_row: dict[str, str] | None,
    expression_row: dict[str, str] | None,
) -> str:
    components = [orthology_component(row)]
    if phase4e_row is not None:
        components.append("PHASE4E_LOCUS_REVIEW")
    if phase5c_row is not None:
        components.append("PHASE5C_REPEAT_CONTEXT_QC")
    if phase6_row is not None:
        components.append("PHASE6_TELOMERE_GENE_READINESS")
    if expression_row is not None:
        components.append("PHASE7E_PHASE8A_RETINA_EXPRESSION_REVIEW")
    return ";".join(components)


def tier1_supported(
    row: dict[str, str],
    artifact_level: str,
    phase4e_row: dict[str, str] | None,
    expression_row: dict[str, str] | None,
) -> bool:
    has_orthology = orthology_component(row) == "ORTHOLOGY_OR_CANDIDATE_MAPPING_PRESENT"
    has_low_artifact = artifact_level in {"low", "not_assessed"}
    has_locus_or_no_duplication_claim = phase4e_row is None or not value_has(phase4e_row.get("copy_number_interpretation_status"), "not_validated")
    expression_ok = expression_row is None or expression_row.get("expression_support_used") == CAUTIOUS_EXPRESSION
    return has_orthology and has_low_artifact and has_locus_or_no_duplication_claim and expression_ok


def final_tier(row: dict[str, str], artifact_level: str, phase4e_row: dict[str, str] | None, expression_row: dict[str, str] | None) -> str:
    current = row["evidence_tier"]
    if current == "Tier 1" and not tier1_supported(row, artifact_level, phase4e_row, expression_row):
        return "Tier 2"
    return current


def reporting_class(final_evidence_tier: str, artifact_level: str, expression_row: dict[str, str] | None) -> str:
    if final_evidence_tier == "Artifact/uncertain":
        return "ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY"
    if final_evidence_tier == "Tier 1" and artifact_level in {"low", "not_assessed"}:
        return "ROBUST_COMPUTATIONAL_EVIDENCE_READY_FOR_REPORTING"
    if final_evidence_tier == "Tier 2" and artifact_level == "high":
        return "PLAUSIBLE_LEAD_ARTIFACT_PRONE"
    if final_evidence_tier == "Tier 2":
        return "PLAUSIBLE_LEAD_REQUIRES_VALIDATION"
    if expression_row and expression_row.get("expression_support_used") == "EXPLORATORY_TISSUE_SIGNAL_ONLY":
        return "EXPLORATORY_SIGNAL_WITH_RETINA_CONTEXT"
    return "EXPLORATORY_SIGNAL"


def audit_status(current_tier: str, final_evidence_tier: str, artifact_level: str, class_label: str) -> str:
    if TIER_ORDER[final_evidence_tier] < TIER_ORDER[current_tier]:
        return "TIER_DOWNGRADED_BY_PHASE8B_RULES"
    if class_label == "ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY":
        return "TIER_RETAINED_AS_ARTIFACT_UNCERTAIN"
    if artifact_level == "high":
        return "TIER_RETAINED_WITH_MAJOR_CAVEATS"
    return "TIER_RETAINED_PENDING_VALIDATION"


def major_caveats(
    row: dict[str, str],
    phase4e_row: dict[str, str] | None,
    phase5c_row: dict[str, str] | None,
    expression_row: dict[str, str] | None,
) -> str:
    caveats: list[str] = []
    artifact = row.get("artifact_risk", "")
    if value_has(artifact, "high"):
        caveats.append("high_artifact_risk")
    if phase4e_row and value_has(phase4e_row.get("copy_number_interpretation_status"), "not_validated"):
        caveats.append("copy_number_not_validated")
    if phase5c_row and value_has(phase5c_row.get("artifact_risk_modifier"), "high_repeat"):
        caveats.append("high_repeat_context_artifact_risk")
    if expression_row and expression_row.get("expression_support_used") != CAUTIOUS_EXPRESSION:
        caveats.append(expression_row.get("expression_support_used", "expression_not_used_as_positive_support"))
    if phase4e_row and value_has(phase4e_row.get("required_validation"), "cross-resource"):
        caveats.append("cross_resource_validation_required")
    if not caveats:
        caveats.append("standard_validation_required")
    return ";".join(dict.fromkeys(caveats))


def required_validation(
    row: dict[str, str],
    phase4e_row: dict[str, str] | None,
    phase5c_row: dict[str, str] | None,
    phase6_row: dict[str, str] | None,
    expression_row: dict[str, str] | None,
) -> str:
    return append_unique_semicolon(
        row.get("required_validation"),
        [
            phase4e_row.get("required_validation", NOT_ASSESSED) if phase4e_row else NOT_ASSESSED,
            phase5c_row.get("required_validation", NOT_ASSESSED) if phase5c_row else NOT_ASSESSED,
            phase6_row.get("required_validation", NOT_ASSESSED) if phase6_row else NOT_ASSESSED,
            expression_row.get("required_validation_update", NOT_ASSESSED) if expression_row else NOT_ASSESSED,
            "Phase 8b final tier audit requires all mechanism-level claims to remain traceable to supporting tables and validation caveats.",
        ],
    )


def supporting_files(
    row: dict[str, str],
    phase4e_row: dict[str, str] | None,
    phase5c_row: dict[str, str] | None,
    phase6_row: dict[str, str] | None,
    expression_row: dict[str, str] | None,
) -> str:
    return append_unique_semicolon(
        row.get("resources_supporting"),
        [
            phase4e_row.get("supporting_files", NOT_ASSESSED) if phase4e_row else NOT_ASSESSED,
            phase5c_row.get("supporting_files", NOT_ASSESSED) if phase5c_row else NOT_ASSESSED,
            phase6_row.get("supporting_files", NOT_ASSESSED) if phase6_row else NOT_ASSESSED,
            expression_row.get("supporting_files", NOT_ASSESSED) if expression_row else NOT_ASSESSED,
            "results/evidence/phase8b_tier_audit.tsv",
        ],
    )


def build_audit_row(
    row: dict[str, str],
    phase4e_row: dict[str, str] | None,
    phase5c_row: dict[str, str] | None,
    phase6_row: dict[str, str] | None,
    expression_row: dict[str, str] | None,
) -> dict[str, str]:
    artifact_level = artifact_risk_level(row, phase5c_row, expression_row)
    final_evidence_tier = final_tier(row, artifact_level, phase4e_row, expression_row)
    class_label = reporting_class(final_evidence_tier, artifact_level, expression_row)
    return {
        "gene_symbol": row["gene_or_pathway"],
        "mechanism": row["mechanism"],
        "current_evidence_tier": row["evidence_tier"],
        "phase8b_final_evidence_tier": final_evidence_tier,
        "tier_change": "UNCHANGED" if final_evidence_tier == row["evidence_tier"] else f"{row['evidence_tier']}->{final_evidence_tier}",
        "phase8b_audit_status": audit_status(row["evidence_tier"], final_evidence_tier, artifact_level, class_label),
        "phase8b_reporting_class": class_label,
        "evidence_components": evidence_components(row, phase4e_row, phase5c_row, phase6_row, expression_row),
        "orthology_component": orthology_component(row),
        "locus_component": locus_component(phase4e_row),
        "repeat_context_component": repeat_component(phase5c_row),
        "telomere_component": telomere_component(phase6_row),
        "expression_component": expression_component(expression_row),
        "artifact_risk_level": artifact_level,
        "major_caveats": major_caveats(row, phase4e_row, phase5c_row, expression_row),
        "required_validation": required_validation(row, phase4e_row, phase5c_row, phase6_row, expression_row),
        "supporting_files": supporting_files(row, phase4e_row, phase5c_row, phase6_row, expression_row),
        "claim_language_guardrail": "Phase 8b is a conservative audit layer. Do not infer mechanism, pathway activity, functional advantage, clinical relevance, causation, validated duplication, telomere length, or organism-wide aging from this table alone.",
    }


def final_interpretation_sentence(audit_row: dict[str, str]) -> str:
    return (
        f"Phase 8b final audit retains {audit_row['phase8b_final_evidence_tier']} and classifies this candidate as "
        f"{audit_row['phase8b_reporting_class']} with {audit_row['artifact_risk_level']} artifact-risk level. "
        "This is a reporting-readiness classification, not a new biological mechanism claim."
    )


def update_evidence_row(row: dict[str, str], audit_row: dict[str, str]) -> dict[str, str]:
    updated = dict(row)
    updated["evidence_tier"] = audit_row["phase8b_final_evidence_tier"]
    updated["evidence_type"] = append_unique_semicolon(row.get("evidence_type"), ["phase8b_final_tier_audit"])
    updated["resources_supporting"] = append_unique_semicolon(row.get("resources_supporting"), [audit_row["supporting_files"]])
    updated["artifact_risk"] = append_unique_semicolon(row.get("artifact_risk"), [f"phase8b_artifact_risk_{audit_row['artifact_risk_level']}"])
    updated["biological_interpretation"] = join_sentences(row.get("biological_interpretation"), final_interpretation_sentence(audit_row))
    updated["relevance_to_aging_longevity"] = join_sentences(
        row.get("relevance_to_aging_longevity"),
        "Phase 8b does not upgrade aging or longevity relevance beyond the validated evidence tier.",
    )
    updated["required_validation"] = append_unique_semicolon(row.get("required_validation"), [audit_row["required_validation"]])
    updated["claim_language_guardrail"] = append_unique_semicolon(row.get("claim_language_guardrail"), [audit_row["claim_language_guardrail"]])
    validate_guardrail_language([updated])
    return updated


def build_mechanism_summary(audit_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in audit_rows:
        grouped[row["mechanism"]].append(row)

    summary_rows: list[dict[str, str]] = []
    for mechanism, rows in sorted(grouped.items()):
        tiers = Counter(row["phase8b_final_evidence_tier"] for row in rows)
        classes = Counter(row["phase8b_reporting_class"] for row in rows)
        expression_supported = sum(1 for row in rows if row["expression_component"] == CAUTIOUS_EXPRESSION)
        repeat_caveats = sum(1 for row in rows if "repeat" in row["major_caveats"].lower())
        locus_review_required = sum(1 for row in rows if value_has(row["major_caveats"], "copy_number_not_validated", "cross_resource"))
        summary_rows.append(
            {
                "mechanism": mechanism,
                "candidate_count": str(len(rows)),
                "tier1_count": str(tiers.get("Tier 1", 0)),
                "tier2_count": str(tiers.get("Tier 2", 0)),
                "tier3_count": str(tiers.get("Tier 3", 0)),
                "artifact_uncertain_count": str(tiers.get("Artifact/uncertain", 0)),
                "robust_ready_count": str(classes.get("ROBUST_COMPUTATIONAL_EVIDENCE_READY_FOR_REPORTING", 0)),
                "plausible_lead_count": str(
                    classes.get("PLAUSIBLE_LEAD_REQUIRES_VALIDATION", 0)
                    + classes.get("PLAUSIBLE_LEAD_ARTIFACT_PRONE", 0)
                ),
                "exploratory_count": str(
                    classes.get("EXPLORATORY_SIGNAL", 0) + classes.get("EXPLORATORY_SIGNAL_WITH_RETINA_CONTEXT", 0)
                ),
                "artifact_prone_count": str(
                    classes.get("PLAUSIBLE_LEAD_ARTIFACT_PRONE", 0)
                    + classes.get("ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY", 0)
                ),
                "expression_supported_count": str(expression_supported),
                "repeat_context_caveat_count": str(repeat_caveats),
                "locus_review_required_count": str(locus_review_required),
                "required_validation": "Use row-level Phase 8b audit entries before making mechanism-level interpretations.",
            }
        )
    return summary_rows


def integrate_phase8b(
    evidence_rows: list[dict[str, str]],
    expression_audit_rows: list[dict[str, str]],
    phase4e_rows: list[dict[str, str]],
    phase5c_rows: list[dict[str, str]],
    phase6_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    if not evidence_rows:
        raise ValueError("Phase 8a integrated evidence has no rows")
    expression_by_gene = index_by(expression_audit_rows, "gene_symbol")
    phase4e_by_gene = index_by(phase4e_rows, "gene_symbol")
    phase5c_by_gene = index_by(phase5c_rows, "gene_symbol")
    phase6_by_gene = index_by(phase6_rows, "gene_symbol")

    audit_rows: list[dict[str, str]] = []
    final_rows: list[dict[str, str]] = []
    for row in evidence_rows:
        gene = row["gene_or_pathway"]
        audit_row = build_audit_row(
            row,
            phase4e_by_gene.get(gene),
            phase5c_by_gene.get(gene),
            phase6_by_gene.get(gene),
            expression_by_gene.get(gene),
        )
        audit_rows.append(audit_row)
        final_rows.append(update_evidence_row(row, audit_row))

    validate_guardrail_language(final_rows)
    return final_rows, audit_rows, build_mechanism_summary(audit_rows)


def write_report(
    final_rows: list[dict[str, str]],
    audit_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    output: Path,
) -> None:
    tiers = Counter(row["phase8b_final_evidence_tier"] for row in audit_rows)
    classes = Counter(row["phase8b_reporting_class"] for row in audit_rows)
    changes = Counter(row["tier_change"] for row in audit_rows)
    lines = [
        "# Phase 8b Final Evidence Scoring And Tier Audit",
        "",
        "Phase 8b consolidates the current evidence streams into a final pre-reporting evidence table and an explicit tier audit.",
        "",
        "## Method Rationale",
        "",
        "The method is deterministic rule-based auditing over validated TSV outputs. This is preferred over a numeric weighted score or machine-learning classifier because the current evidence streams are categorical, heterogeneous, and incomplete. A rule-based audit keeps each retained or adjusted tier traceable to orthology, locus validation, repeat context, telomere readiness, and retina-expression caveats.",
        "",
        "## Final Tier Summary",
        "",
        f"- Evidence rows audited: {len(audit_rows)}",
        f"- Final tier counts: {dict(sorted(tiers.items()))}",
        f"- Reporting-class counts: {dict(sorted(classes.items()))}",
        f"- Tier-change counts: {dict(sorted(changes.items()))}",
        "",
        "## Mechanism Summary",
        "",
        "| Mechanism | Candidates | Tier 2 | Artifact/uncertain | Expression-supported | Repeat caveat | Locus review required |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["mechanism"],
                    row["candidate_count"],
                    row["tier2_count"],
                    row["artifact_uncertain_count"],
                    row["expression_supported_count"],
                    row["repeat_context_caveat_count"],
                    row["locus_review_required_count"],
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Guardrail",
            "",
            "Phase 8b is a pre-reporting audit. It does not create mechanism claims, pathway-state claims, telomere-length claims, validated duplication claims, causation claims, or human translational claims.",
            "",
            "Supporting tables:",
            "",
            "- `results/evidence/phase8b_final_integrated_evidence.tsv`",
            "- `results/evidence/phase8b_tier_audit.tsv`",
            "- `results/evidence/phase8b_mechanism_summary.tsv`",
        ]
    )
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase8b(
    config_path: Path,
    final_evidence_output: Path,
    tier_audit_output: Path,
    mechanism_summary_output: Path,
    report_output: Path,
) -> None:
    config = read_yaml(config_path)
    if "phase8b_final_evidence_scoring" not in config:
        raise ValueError("config.yaml is missing phase8b_final_evidence_scoring")
    phase = config["phase8b_final_evidence_scoring"]
    evidence_rows = read_tsv(Path(phase["phase8a_integrated_evidence"]), EVIDENCE_COLUMNS)
    expression_audit_rows = read_tsv(Path(phase["phase8a_expression_integration_audit"]), PHASE8A_INTEGRATION_AUDIT_COLUMNS)
    phase4e_rows = read_tsv(Path(phase["phase4e_gene_hardened_summary"]), PHASE4E_GENE_HARDENING_COLUMNS)
    phase5c_rows = read_tsv(Path(phase["phase5c_gene_repeat_qc_summary"]), PHASE5C_GENE_QC_COLUMNS)
    phase6_rows = read_tsv(Path(phase["phase6_telomere_gene_audit"]), PHASE6_TELOMERE_GENE_AUDIT_COLUMNS)

    final_rows, audit_rows, summary_rows = integrate_phase8b(
        evidence_rows,
        expression_audit_rows,
        phase4e_rows,
        phase5c_rows,
        phase6_rows,
    )
    write_tsv(final_evidence_output, final_rows, EVIDENCE_COLUMNS)
    write_tsv(tier_audit_output, audit_rows, PHASE8B_TIER_AUDIT_COLUMNS)
    write_tsv(mechanism_summary_output, summary_rows, PHASE8B_MECHANISM_SUMMARY_COLUMNS)
    write_report(final_rows, audit_rows, summary_rows, report_output)
    LOGGER.info("Wrote Phase 8b final evidence audit for %d evidence rows", len(audit_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Phase 8b final integrated evidence and tier audit.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--final-evidence-output", type=Path, required=True)
    parser.add_argument("--tier-audit-output", type=Path, required=True)
    parser.add_argument("--mechanism-summary-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_phase8b(
        args.config,
        args.final_evidence_output,
        args.tier_audit_output,
        args.mechanism_summary_output,
        args.report_output,
    )


if __name__ == "__main__":
    main()
