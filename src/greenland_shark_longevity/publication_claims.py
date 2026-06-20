"""Register and audit published claims against current local evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

from .schemas import ALLOWED_MECHANISMS, PUBLICATION_CLAIM_AUDIT_COLUMNS, REFERENCE_FILE_COLUMNS
from .utils import read_tsv, read_yaml, write_tsv


REQUIRED_CLAIM_FIELDS = {
    "claim_id",
    "publication_id",
    "doi",
    "resource_id",
    "mechanism",
    "gene_or_feature",
    "reported_signal",
    "required_evidence",
    "current_repo_status",
    "interpretation_guardrail",
    "next_validation",
}

EVIDENCE_TO_FILE_ROLES = {
    "assembly_stats": {"assembly_stats"},
    "assembly_report": {"assembly_report"},
    "genome_fasta_checksum": {"md5checksums"},
    "genome_wide_protein_set": {"protein_fasta"},
    "annotation_coordinates": {"annotation_gff", "annotation_genbank"},
    "raw_reads_or_variant_calls": {"sra_reads", "variant_calls"},
}


def _validate_claims(claim_config: dict) -> list[dict]:
    claims = claim_config.get("claims")
    if not isinstance(claims, list) or not claims:
        raise ValueError("Publication claim YAML must contain a non-empty claims list")
    seen: set[str] = set()
    for claim in claims:
        if not isinstance(claim, dict):
            raise ValueError("Each publication claim must be a mapping")
        missing = sorted(REQUIRED_CLAIM_FIELDS - set(claim))
        if missing:
            raise ValueError(f"Claim {claim.get('claim_id', '<unknown>')} missing fields: {', '.join(missing)}")
        if claim["claim_id"] in seen:
            raise ValueError(f"Duplicate publication claim_id: {claim['claim_id']}")
        seen.add(claim["claim_id"])
        if claim["mechanism"] not in ALLOWED_MECHANISMS:
            raise ValueError(f"Claim {claim['claim_id']} has unsupported mechanism {claim['mechanism']!r}")
        if not isinstance(claim["required_evidence"], list) or not claim["required_evidence"]:
            raise ValueError(f"Claim {claim['claim_id']} must include required_evidence list")
    return claims


def _inventory_by_resource(inventory_rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in inventory_rows:
        grouped.setdefault(row["resource_id"], []).append(row)
    return grouped


def _available_roles(rows: list[dict[str, str]]) -> set[str]:
    return {row["file_role"] for row in rows if row["file_status"] == "LOCAL_AVAILABLE"}


def _registered_roles(rows: list[dict[str, str]]) -> set[str]:
    unavailable = {"NOT_AVAILABLE", "NOT_AVAILABLE_FROM_REGISTERED_SOURCE"}
    roles = set()
    for row in rows:
        if row["url"] not in unavailable and row["file_status"] != "NOT_AVAILABLE_FROM_REGISTERED_SOURCE":
            roles.add(row["file_role"])
    return roles


def _blockers(required_evidence: list[str], rows: list[dict[str, str]]) -> list[str]:
    available = _available_roles(rows)
    registered = _registered_roles(rows)
    blockers = []
    for evidence in required_evidence:
        roles = EVIDENCE_TO_FILE_ROLES.get(evidence)
        if roles is None:
            blockers.append(f"{evidence}:REQUIRES_NEW_WORKFLOW_OR_MANUAL_REVIEW")
            continue
        if roles & available:
            continue
        if roles & registered:
            blockers.append(f"{evidence}:REGISTERED_NOT_DOWNLOADED_OR_NOT_PARSED")
        else:
            blockers.append(f"{evidence}:NOT_AVAILABLE_FROM_CURRENT_REGISTERED_RESOURCE")
    return blockers


def _tier_for_claim(claim: dict, blockers: list[str]) -> str:
    status = str(claim["current_repo_status"])
    if status == "PARTIALLY_REPRODUCED_RESOURCE_QC":
        return "Resource-quality observation"
    if blockers:
        return "NOT_ASSESSED"
    return "Tier 2"


def audit_publication_claims(claims_path: Path, inventory_path: Path) -> list[dict[str, str]]:
    claims = _validate_claims(read_yaml(claims_path))
    inventory = _inventory_by_resource(read_tsv(inventory_path, REFERENCE_FILE_COLUMNS))
    rows: list[dict[str, str]] = []
    for claim in claims:
        resource_rows = inventory.get(claim["resource_id"], [])
        blockers = _blockers([str(item) for item in claim["required_evidence"]], resource_rows)
        available = sorted(_available_roles(resource_rows))
        rows.append(
            {
                "claim_id": str(claim["claim_id"]),
                "publication_id": str(claim["publication_id"]),
                "doi": str(claim["doi"]),
                "resource_id": str(claim["resource_id"]),
                "mechanism": str(claim["mechanism"]),
                "gene_or_feature": str(claim["gene_or_feature"]),
                "reported_signal": str(claim["reported_signal"]),
                "current_repo_status": str(claim["current_repo_status"]),
                "evidence_tier": _tier_for_claim(claim, blockers),
                "available_supporting_files": ";".join(available) if available else "NONE",
                "blockers": ";".join(blockers) if blockers else "NONE",
                "next_validation": str(claim["next_validation"]),
                "interpretation_guardrail": str(claim["interpretation_guardrail"]),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit registered publication claims against current local resources.")
    parser.add_argument("--claims", type=Path, required=True)
    parser.add_argument("--reference-inventory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = audit_publication_claims(args.claims, args.reference_inventory)
    write_tsv(args.output, rows, PUBLICATION_CLAIM_AUDIT_COLUMNS)


if __name__ == "__main__":
    main()
