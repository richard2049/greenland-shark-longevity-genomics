"""Evidence-tier scoring."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .schemas import DUPLICATION_AUDIT_COLUMNS, EVIDENCE_COLUMNS, GUARDRAIL_TERMS
from .utils import as_bool, read_tsv, write_tsv


def score_audit_row(row: dict[str, str]) -> dict[str, str]:
    demo_only = as_bool(row.get("demo_only"))
    orthology = row.get("orthology_support") == "yes"
    complete_domain = row.get("domain_integrity") == "complete"
    separable_loci = row.get("separable_loci") == "yes"
    coordinate_support = row.get("coordinate_support") == "yes"
    cross_resource = row.get("cross_resource_support") == "yes"
    artifact_risk = row.get("artifact_risk", "not_assessable")
    copy_count = int(row.get("copy_count", "0") or 0)

    if demo_only:
        tier = "Tier 3"
        interpretation = "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE; schema example only."
        aging_relevance = "NOT_ASSESSED"
    elif artifact_risk == "high" or not orthology:
        tier = "Artifact/uncertain"
        interpretation = "Do not claim biologically; evidence is incomplete or artifact-prone."
        aging_relevance = "NOT_ASSESSED"
    elif orthology and complete_domain and (cross_resource or (separable_loci and coordinate_support)):
        tier = "Tier 1"
        interpretation = "Robust computational evidence under the current schema."
        aging_relevance = "Potentially relevant only within the stated mechanism and validation limits."
    elif orthology and (complete_domain or coordinate_support):
        tier = "Tier 2"
        interpretation = "Plausible but incomplete computational evidence."
        aging_relevance = "Hypothesis-generating; requires validation before biological interpretation."
    else:
        tier = "Tier 3"
        interpretation = "Exploratory signal only."
        aging_relevance = "NOT_ASSESSED"

    if tier == "Tier 1" and not (orthology and complete_domain and (cross_resource or (separable_loci and coordinate_support))):
        raise ValueError(f"Invalid Tier 1 assignment for {row.get('gene_symbol')}")

    evidence_type = "candidate_copy_number" if copy_count > 1 else "candidate_presence_or_not_assessed"
    guardrail = (
        "Avoid activated/inactivated/absent/causal wording; demo rows are not biological evidence."
        if demo_only
        else "Avoid activated/inactivated/absent/causal wording unless direct evidence supports it."
    )
    return {
        "mechanism": row["mechanism"],
        "gene_or_pathway": row["gene_symbol"],
        "evidence_type": evidence_type,
        "evidence_tier": tier,
        "resources_supporting": row.get("resources_supporting", "NOT_ASSESSED"),
        "artifact_risk": artifact_risk,
        "biological_interpretation": interpretation,
        "relevance_to_aging_longevity": aging_relevance,
        "translational_relevance": "NOT_ASSESSED",
        "required_validation": row.get("required_validation", "REQUIRES_VALIDATION"),
        "claim_language_guardrail": guardrail,
    }


def validate_guardrail_language(rows: list[dict[str, str]]) -> None:
    checked_columns = ["biological_interpretation", "relevance_to_aging_longevity", "translational_relevance"]
    for row in rows:
        text = " ".join(row.get(column, "") for column in checked_columns).lower()
        for term in GUARDRAIL_TERMS:
            if term in text and "avoid" not in text:
                raise ValueError(f"Guardrail term {term!r} appears in evidence row for {row.get('gene_or_pathway')}")


def score_evidence(audit_path: Path) -> list[dict[str, str]]:
    audit_rows = read_tsv(audit_path, DUPLICATION_AUDIT_COLUMNS)
    evidence_rows = [score_audit_row(row) for row in audit_rows]
    validate_guardrail_language(evidence_rows)
    return evidence_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Score integrated evidence tiers.")
    parser.add_argument("--duplication-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    write_tsv(args.output, score_evidence(args.duplication_audit), EVIDENCE_COLUMNS)


if __name__ == "__main__":
    main()

