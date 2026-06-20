"""Phase 8a integration of hardened retina expression support into evidence tables."""

from __future__ import annotations

import argparse
import logging
from collections import Counter
from pathlib import Path

from .evidence import validate_guardrail_language
from .schemas import EVIDENCE_COLUMNS, PHASE7E_CANDIDATE_HARDENED_COLUMNS, PHASE8A_INTEGRATION_AUDIT_COLUMNS
from .utils import NOT_ASSESSED, ensure_parent, read_tsv, read_yaml, write_tsv

LOGGER = logging.getLogger(__name__)

CAUTIOUS_SUPPORT = "CAUTIOUS_RETINA_EXPRESSION_SUPPORT_FOR_PHASE8"
LIMITED_SUPPORT = "LIMITED_EXPLORATORY_EXPRESSION_SUPPORT"
DEFER_SUPPORT = "DEFER_PHASE8_EXPRESSION_SUPPORT_UNTIL_LOCUS_REVIEW"
DO_NOT_USE_SUPPORT = "DO_NOT_USE_AS_PHASE8_EXPRESSION_SUPPORT"
NO_SUPPORT = "NO_RETINA_EXPRESSION_SUPPORT_UNDER_CURRENT_REFERENCE"


def append_unique_semicolon(value: str | None, additions: list[str]) -> str:
    parts: list[str] = []
    for item in [value or "", *additions]:
        if item in {"", NOT_ASSESSED, None}:
            continue
        for part in str(item).split(";"):
            cleaned = part.strip()
            if cleaned and cleaned != NOT_ASSESSED and cleaned not in parts:
                parts.append(cleaned)
    return ";".join(parts) if parts else NOT_ASSESSED


def join_sentences(*sentences: str | None) -> str:
    cleaned = [str(sentence).strip() for sentence in sentences if sentence and str(sentence).strip() != NOT_ASSESSED]
    if not cleaned:
        return NOT_ASSESSED
    return " ".join(sentence if sentence.endswith(".") else f"{sentence}." for sentence in cleaned)


def expression_evidence_tier(row: dict[str, str]) -> str:
    support = row["phase7e_expression_support_level"]
    if support == CAUTIOUS_SUPPORT:
        return "Tier 2"
    if support == LIMITED_SUPPORT:
        return "Tier 3"
    return "Artifact/uncertain"


def expression_evidence_type(row: dict[str, str]) -> str:
    support = row["phase7e_expression_support_level"]
    if support == CAUTIOUS_SUPPORT:
        return "phase8a_retina_expression_support_cautious"
    if support == LIMITED_SUPPORT:
        return "phase8a_retina_expression_support_limited_exploratory"
    if support == DEFER_SUPPORT:
        return "phase8a_expression_support_deferred_locus_review"
    if support == NO_SUPPORT:
        return "phase8a_no_retina_expression_support_current_reference"
    if support == DO_NOT_USE_SUPPORT:
        return "phase8a_expression_not_interpretable"
    return "phase8a_expression_not_assessed"


def expression_support_used(row: dict[str, str]) -> str:
    support = row["phase7e_expression_support_level"]
    if support == CAUTIOUS_SUPPORT:
        return "CAUTIOUS_TISSUE_SPECIFIC_SUPPORT"
    if support == LIMITED_SUPPORT:
        return "EXPLORATORY_TISSUE_SIGNAL_ONLY"
    return "NO_POSITIVE_EXPRESSION_SUPPORT_USED"


def integration_action(row: dict[str, str]) -> str:
    support = row["phase7e_expression_support_level"]
    if support == CAUTIOUS_SUPPORT:
        return "APPEND_CAUTION_RETINA_EXPRESSION_SUPPORT_NO_TIER_CHANGE"
    if support == LIMITED_SUPPORT:
        return "APPEND_LIMITED_EXPLORATORY_EXPRESSION_CAVEAT_NO_TIER_CHANGE"
    if support == DEFER_SUPPORT:
        return "APPEND_DEFERRED_EXPRESSION_CAVEAT_NO_TIER_CHANGE"
    if support == NO_SUPPORT:
        return "APPEND_NO_RETINA_EXPRESSION_SUPPORT_CAVEAT_NO_TIER_CHANGE"
    return "APPEND_NOT_INTERPRETABLE_EXPRESSION_CAVEAT_NO_TIER_CHANGE"


def artifact_risk_update(row: dict[str, str]) -> str:
    support = row["phase7e_expression_support_level"]
    support_label = {
        CAUTIOUS_SUPPORT: "phase7e_retina_expression_qc_caveat",
        LIMITED_SUPPORT: "limited_exploratory_retina_expression_signal",
        DEFER_SUPPORT: "expression_support_deferred_locus_or_reference_ambiguity",
        NO_SUPPORT: "no_retina_expression_support_current_reference",
        DO_NOT_USE_SUPPORT: "expression_not_interpretable_current_reference",
    }.get(support, "expression_support_not_assessed")
    return append_unique_semicolon(row.get("artifact_risk"), [support_label])


def required_validation_update(row: dict[str, str]) -> str:
    return append_unique_semicolon(
        row.get("required_validation"),
        [
            "Do not upgrade biological evidence tier from expression alone.",
            "Require orthology, locus identity, domain integrity, cross-resource support, and expression-reference validation before stronger claims.",
        ],
    )


def expression_interpretation(row: dict[str, str]) -> str:
    gene = row["gene_symbol"]
    support = row["phase7e_expression_support_level"]
    detected = row.get("detected_run_count", NOT_ASSESSED)
    quantified = row.get("quantified_run_count", NOT_ASSESSED)
    median_tpm = row.get("median_tpm", NOT_ASSESSED)
    if support == CAUTIOUS_SUPPORT:
        return (
            f"Phase 7e records cautious retina-specific expression detection for {gene} in "
            f"{detected}/{quantified} quantified runs (median TPM {median_tpm}). "
            "This supports transcript detection in retina only and does not establish differential expression, pathway state, function, or longevity mechanism."
        )
    if support == LIMITED_SUPPORT:
        return (
            f"Phase 7e records limited retina expression signal for {gene} in {detected}/{quantified} quantified runs "
            f"(median TPM {median_tpm}). Treat this as exploratory tissue signal only."
        )
    if support == DEFER_SUPPORT:
        return (
            f"Phase 7e does not use the quantified retina signal for {gene} as expression support because "
            "locus or reference ambiguity requires review."
        )
    if support == NO_SUPPORT:
        return (
            f"Phase 7e does not provide retina expression support for {gene} under the current thresholds and reference. "
            "This is condition-specific and should not be treated as gene loss or lack of function."
        )
    return f"Phase 7e marks {gene} expression as not interpretable under the current reference and validation state."


def build_expression_evidence_row(row: dict[str, str]) -> dict[str, str]:
    evidence_row = {
        "mechanism": row["mechanism"],
        "gene_or_pathway": row["gene_symbol"],
        "evidence_type": expression_evidence_type(row),
        "evidence_tier": expression_evidence_tier(row),
        "resources_supporting": append_unique_semicolon(
            row.get("supporting_files"),
            [
                "results/rnaseq/phase7e_candidate_expression_hardened.tsv",
                "results/rnaseq/phase7e_run_qc_review.tsv",
                "results/rnaseq/phase7e_parameter_review.tsv",
            ],
        ),
        "artifact_risk": artifact_risk_update(row),
        "biological_interpretation": expression_interpretation(row),
        "relevance_to_aging_longevity": "Retina-specific transcript detection can support tissue-expression availability only; it is not aging or longevity evidence by itself.",
        "translational_relevance": "NOT_ASSESSED",
        "required_validation": required_validation_update(row),
        "claim_language_guardrail": "Expression support is retina-specific and cannot upgrade biological tiers without independent validation; do not infer activation, differential expression, pathway state, organism-wide aging, causation, functional advantage, gene loss, or longevity mechanism.",
    }
    validate_guardrail_language([evidence_row])
    return evidence_row


def update_integrated_evidence_row(base_row: dict[str, str], expression_row: dict[str, str]) -> dict[str, str]:
    updated = dict(base_row)
    updated["evidence_tier"] = base_row["evidence_tier"]
    updated["evidence_type"] = append_unique_semicolon(base_row.get("evidence_type"), [expression_evidence_type(expression_row)])
    updated["resources_supporting"] = append_unique_semicolon(
        base_row.get("resources_supporting"),
        [
            expression_row.get("supporting_files", NOT_ASSESSED),
            "results/rnaseq/phase7e_candidate_expression_hardened.tsv",
            "results/rnaseq/phase7e_run_qc_review.tsv",
            "results/rnaseq/phase7e_parameter_review.tsv",
        ],
    )
    updated["artifact_risk"] = append_unique_semicolon(base_row.get("artifact_risk"), [artifact_risk_update(expression_row)])
    updated["biological_interpretation"] = join_sentences(base_row.get("biological_interpretation"), expression_interpretation(expression_row))
    updated["relevance_to_aging_longevity"] = join_sentences(
        base_row.get("relevance_to_aging_longevity"),
        "Phase 7e adds retina-specific expression context only; the integrated biological tier is unchanged.",
    )
    updated["required_validation"] = append_unique_semicolon(base_row.get("required_validation"), [required_validation_update(expression_row)])
    updated["claim_language_guardrail"] = append_unique_semicolon(
        base_row.get("claim_language_guardrail"),
        [
            expression_row.get("claim_language_guardrail", NOT_ASSESSED),
            "Phase 8a expression context cannot upgrade biological tiers without all required validation criteria.",
        ],
    )
    if updated["evidence_tier"] != base_row["evidence_tier"]:
        raise ValueError(f"Phase 8a attempted to change evidence tier for {base_row['gene_or_pathway']}")
    validate_guardrail_language([updated])
    return updated


def build_audit_row(base_row: dict[str, str] | None, expression_row: dict[str, str], integrated_tier: str) -> dict[str, str]:
    base_tier = base_row["evidence_tier"] if base_row else NOT_ASSESSED
    return {
        "gene_symbol": expression_row["gene_symbol"],
        "base_evidence_tier": base_tier,
        "expression_evidence_tier": expression_evidence_tier(expression_row),
        "integrated_evidence_tier": integrated_tier,
        "tier_change": "UNCHANGED" if base_row and base_tier == integrated_tier else "NO_BASE_ROW_ADDED_AS_UNCERTAIN",
        "phase7e_expression_support_level": expression_row["phase7e_expression_support_level"],
        "expression_interpretation_status": expression_row["expression_interpretation_status"],
        "integration_action": integration_action(expression_row),
        "expression_support_used": expression_support_used(expression_row),
        "artifact_risk_update": artifact_risk_update(expression_row),
        "required_validation_update": required_validation_update(expression_row),
        "supporting_files": append_unique_semicolon(
            expression_row.get("supporting_files"),
            ["results/evidence/phase8a_expression_support_evidence.tsv"],
        ),
    }


def integrate_phase8a_expression_support(
    base_rows: list[dict[str, str]],
    expression_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    if not base_rows:
        raise ValueError("Base evidence table has no rows")
    if not expression_rows:
        raise ValueError("Phase 7e candidate expression table has no rows")

    expression_by_gene = {row["gene_symbol"]: row for row in expression_rows}
    expression_evidence_rows = [build_expression_evidence_row(row) for row in sorted(expression_rows, key=lambda item: item["gene_symbol"])]
    integrated_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    integrated_genes: set[str] = set()

    for base_row in base_rows:
        gene = base_row["gene_or_pathway"]
        expression_row = expression_by_gene.get(gene)
        if expression_row is None:
            integrated_rows.append(base_row)
            continue
        updated = update_integrated_evidence_row(base_row, expression_row)
        integrated_rows.append(updated)
        audit_rows.append(build_audit_row(base_row, expression_row, updated["evidence_tier"]))
        integrated_genes.add(gene)

    base_genes = {row["gene_or_pathway"] for row in base_rows}
    for expression_row in sorted(expression_rows, key=lambda item: item["gene_symbol"]):
        gene = expression_row["gene_symbol"]
        if gene in base_genes:
            continue
        fallback = build_expression_evidence_row(expression_row)
        fallback["evidence_tier"] = "Artifact/uncertain"
        fallback["artifact_risk"] = append_unique_semicolon(fallback["artifact_risk"], ["missing_base_evidence_row"])
        fallback["biological_interpretation"] = join_sentences(
            fallback["biological_interpretation"],
            "Phase 8a cannot integrate this expression row with upstream evidence because no matching base evidence row was found.",
        )
        fallback["required_validation"] = append_unique_semicolon(
            fallback["required_validation"],
            ["Create or validate upstream orthology/candidate evidence before using this expression row."],
        )
        integrated_rows.append(fallback)
        audit_rows.append(build_audit_row(None, expression_row, fallback["evidence_tier"]))
        integrated_genes.add(gene)

    validate_guardrail_language(expression_evidence_rows)
    validate_guardrail_language(integrated_rows)
    return expression_evidence_rows, audit_rows, integrated_rows


def write_report(
    expression_evidence_rows: list[dict[str, str]],
    audit_rows: list[dict[str, str]],
    integrated_rows: list[dict[str, str]],
    output: Path,
) -> None:
    support_counts = Counter(row["phase7e_expression_support_level"] for row in audit_rows)
    expression_tier_counts = Counter(row["evidence_tier"] for row in expression_evidence_rows)
    integrated_tier_counts = Counter(row["evidence_tier"] for row in integrated_rows)
    tier_changes = Counter(row["tier_change"] for row in audit_rows)
    lines = [
        "# Phase 8a Expression Evidence Integration",
        "",
        "Phase 8a integrates hardened Phase 7e retina expression support into the evidence-scoring model without upgrading biological evidence tiers.",
        "",
        "## Method Rationale",
        "",
        "The method is deterministic rule-based table integration. This is preferred over a new statistical score because the current RNA-seq layer is retina-only, has three runs, has no defensible differential-expression contrast, and includes reference/mapping caveats. Expression can add tissue-specific support or caveats, but it cannot upgrade a biological tier unless all upstream validation criteria are met in a later scoring step.",
        "",
        "## Audit Summary",
        "",
        f"- Phase 7e candidate rows integrated or reviewed: {len(audit_rows)}",
        f"- Expression support evidence rows written: {len(expression_evidence_rows)}",
        f"- Integrated evidence rows written: {len(integrated_rows)}",
        f"- Tier-change audit: {dict(sorted(tier_changes.items()))}",
        f"- Expression evidence tiers: {dict(sorted(expression_tier_counts.items()))}",
        f"- Integrated evidence tiers: {dict(sorted(integrated_tier_counts.items()))}",
        "",
        "## Phase 7e Support Levels",
        "",
    ]
    for support, count in sorted(support_counts.items()):
        lines.append(f"- `{support}`: {count}")
    lines.extend(
        [
            "",
            "## Guardrail",
            "",
            "Phase 8a supports cautious retina-specific expression context only. It does not support pathway activity, differential expression, whole-organism aging interpretation, causation, functional advantage, gene loss, or longevity mechanism.",
            "",
            "Supporting tables:",
            "",
            "- `results/rnaseq/phase7e_candidate_expression_hardened.tsv`",
            "- `results/evidence/phase8a_expression_support_evidence.tsv`",
            "- `results/evidence/phase8a_expression_integration_audit.tsv`",
            "- `results/evidence/phase8a_integrated_evidence.tsv`",
        ]
    )
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase8a(
    config_path: Path,
    expression_evidence_output: Path,
    integration_audit_output: Path,
    integrated_output: Path,
    report_output: Path,
) -> None:
    config = read_yaml(config_path)
    if "phase8a_expression_integration" not in config:
        raise ValueError("config.yaml is missing phase8a_expression_integration")
    phase = config["phase8a_expression_integration"]
    base_rows = read_tsv(Path(phase["base_evidence"]), EVIDENCE_COLUMNS)
    expression_rows = read_tsv(Path(phase["phase7e_candidate_expression"]), PHASE7E_CANDIDATE_HARDENED_COLUMNS)
    expression_evidence_rows, audit_rows, integrated_rows = integrate_phase8a_expression_support(base_rows, expression_rows)
    write_tsv(expression_evidence_output, expression_evidence_rows, EVIDENCE_COLUMNS)
    write_tsv(integration_audit_output, audit_rows, PHASE8A_INTEGRATION_AUDIT_COLUMNS)
    write_tsv(integrated_output, integrated_rows, EVIDENCE_COLUMNS)
    write_report(expression_evidence_rows, audit_rows, integrated_rows, report_output)
    LOGGER.info("Wrote Phase 8a expression evidence integration for %d candidate rows", len(audit_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Integrate Phase 7e expression support into Phase 8a evidence tables.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--expression-evidence-output", type=Path, required=True)
    parser.add_argument("--integration-audit-output", type=Path, required=True)
    parser.add_argument("--integrated-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_phase8a(args.config, args.expression_evidence_output, args.integration_audit_output, args.integrated_output, args.report_output)


if __name__ == "__main__":
    main()
