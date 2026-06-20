"""Phase 3 public-input discovery and OrthoFinder readiness checks."""

from __future__ import annotations

import argparse
import gzip
import logging
import shutil
from pathlib import Path

from .schemas import (
    ORTHOFINDER_INPUT_MANIFEST_COLUMNS,
    PHASE3_SOURCE_CHECK_COLUMNS,
    REFERENCE_FILE_COLUMNS,
    REFERENCE_GENE_COORDINATE_COLUMNS,
    REFERENCE_PROTEIN_QC_COLUMNS,
)
from .utils import as_bool, configure_logging, read_tsv, read_yaml, write_tsv

NOT_AVAILABLE_VALUES = {"", "TODO", "NOT_ASSESSED", "NOT_AVAILABLE"}
SOURCE_CHECK_REQUIRED_FIELDS = {
    "source_id",
    "resource_id",
    "source_type",
    "url_or_accession",
    "searched_for",
    "status",
    "found_file_types",
    "selected_for_download",
    "local_path",
    "notes",
}


def _as_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        return ";".join(str(item) for item in value)
    return str(value)


def configured_source_checks(config: dict) -> list[dict[str, str]]:
    retrieval_date = _as_text(config.get("phase3_input_discovery", {}).get("retrieval_date", "NOT_ASSESSED"))
    source_checks = config.get("phase3_input_discovery", {}).get("source_checks", [])
    if not isinstance(source_checks, list) or not source_checks:
        raise ValueError("Config must include phase3_input_discovery.source_checks")

    rows = []
    seen: set[str] = set()
    for check in source_checks:
        if not isinstance(check, dict):
            raise ValueError("Each phase3 source check must be a mapping")
        missing = sorted(SOURCE_CHECK_REQUIRED_FIELDS - set(check))
        if missing:
            raise ValueError(f"Source check {check.get('source_id', '<unknown>')} missing fields: {', '.join(missing)}")
        source_id = str(check["source_id"])
        if source_id in seen:
            raise ValueError(f"Duplicate phase3 source_id: {source_id}")
        seen.add(source_id)
        row = {column: _as_text(check.get(column, "")) for column in PHASE3_SOURCE_CHECK_COLUMNS}
        row["retrieval_date"] = row["retrieval_date"] or retrieval_date
        row["selected_for_download"] = str(as_bool(check.get("selected_for_download")))
        rows.append(row)
    return rows


def _species_by_resource(config: dict) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for species in config.get("species", []):
        species_id = str(species.get("species_id", "NOT_ASSESSED"))
        for resource_id in species.get("resource_ids", []):
            mapping[str(resource_id)] = species_id
    return mapping


def _resource_by_id(config: dict) -> dict[str, dict]:
    return {str(resource["resource_id"]): resource for resource in config.get("resources", [])}


def _first_by_role(rows: list[dict[str, str]], resource_id: str, role: str) -> dict[str, str] | None:
    for row in rows:
        if row["resource_id"] == resource_id and row["file_role"] == role:
            return row
    return None


def _protein_qc_by_resource(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["resource_id"]: row for row in rows}


def _is_missing(value: str) -> bool:
    return value in NOT_AVAILABLE_VALUES


def _protein_scope(resource: dict, protein_row: dict[str, str] | None) -> str:
    if protein_row is None or _is_missing(protein_row.get("url", "")):
        return "NO_PROTEIN_FASTA"
    notes = protein_row.get("notes", "").lower()
    resource_type = str(resource.get("resource_type", "")).lower()
    if "not a whole-genome protein set" in notes or "coding-sequence" in notes or "cds-derived" in notes:
        return "CDS_SUBSET_NOT_GENOME_WIDE"
    if "raw_reads_and_coding_sequences" in resource_type:
        return "CDS_SUBSET_NOT_GENOME_WIDE"
    return "GENOME_WIDE_PROTEIN_FASTA_CANDIDATE"


def _orthofinder_blocker(protein_row: dict[str, str] | None, scope: str, qc_row: dict[str, str] | None) -> str:
    if scope == "NO_PROTEIN_FASTA":
        return "NO_GENOME_WIDE_PROTEIN_FASTA_REGISTERED"
    if scope != "GENOME_WIDE_PROTEIN_FASTA_CANDIDATE":
        return "PROTEIN_FASTA_IS_NOT_GENOME_WIDE"
    if protein_row is None or protein_row.get("file_status") != "LOCAL_AVAILABLE":
        return "PROTEIN_FASTA_REGISTERED_BUT_NOT_LOCAL"
    if qc_row is None or qc_row.get("status") != "LOCAL_AVAILABLE":
        return "PROTEIN_QC_NOT_AVAILABLE"
    return "NONE"


def _safe_orthofinder_filename(species_id: str, resource_id: str) -> str:
    cleaned = f"{species_id}__{resource_id}.faa"
    return "".join(char if char.isalnum() or char in {"_", "-", "."} else "_" for char in cleaned)


def stage_orthofinder_input(source_path: Path, target_path: Path) -> None:
    target_path.parent.mkdir(parents=True, exist_ok=True)
    if source_path.suffix == ".gz":
        with gzip.open(source_path, "rb") as source, target_path.open("wb") as target:
            shutil.copyfileobj(source, target)
        return
    shutil.copy2(source_path, target_path)


def build_orthofinder_input_manifest(
    config: dict,
    inventory_rows: list[dict[str, str]],
    protein_qc_rows: list[dict[str, str]],
    copy_inputs: bool = False,
) -> list[dict[str, str]]:
    species_by_resource = _species_by_resource(config)
    resources = _resource_by_id(config)
    qc_by_resource = _protein_qc_by_resource(protein_qc_rows)
    input_dir = Path(config.get("orthofinder", {}).get("input_dir", "data/interim/orthofinder_input"))
    rows: list[dict[str, str]] = []

    for resource_id in sorted(resources):
        resource = resources[resource_id]
        protein_row = _first_by_role(inventory_rows, resource_id, "protein_fasta")
        qc_row = qc_by_resource.get(resource_id)
        scope = _protein_scope(resource, protein_row)
        blocker = _orthofinder_blocker(protein_row, scope, qc_row)
        species_id = species_by_resource.get(resource_id, "NOT_ASSESSED")
        target_path = "NOT_PREPARED"
        ready = blocker == "NONE"
        source_path = protein_row["local_path"] if protein_row else "NOT_AVAILABLE"
        protein_count = qc_row.get("protein_count", "NOT_ASSESSED") if qc_row else "NOT_ASSESSED"
        qc_status = qc_row.get("status", "NOT_ASSESSED") if qc_row else "NOT_ASSESSED"

        if ready:
            target = input_dir / _safe_orthofinder_filename(species_id, resource_id)
            target_path = str(target)
            if copy_inputs:
                stage_orthofinder_input(Path(source_path), target)
                logging.info("Prepared OrthoFinder input %s", target)

        rows.append(
            {
                "species_id": species_id,
                "resource_id": resource_id,
                "source_protein_fasta": source_path if not _is_missing(source_path) else "NOT_AVAILABLE",
                "orthofinder_input_path": target_path,
                "protein_set_scope": scope,
                "protein_count": protein_count,
                "qc_status": qc_status,
                "orthofinder_ready": str(ready),
                "blocker": blocker,
                "next_action": _next_action_for_blocker(blocker),
                "notes": _orthofinder_notes(scope, protein_row),
            }
        )
    return rows


def _next_action_for_blocker(blocker: str) -> str:
    actions = {
        "NONE": "Ready as OrthoFinder input; execution remains disabled until the run rule is explicitly added/enabled.",
        "NO_GENOME_WIDE_PROTEIN_FASTA_REGISTERED": "Locate an author-provided protein set or generate annotation/rescue outputs in a controlled future workflow.",
        "PROTEIN_FASTA_IS_NOT_GENOME_WIDE": "Keep this FASTA out of OrthoFinder; use only for sequence-level provenance/QC.",
        "PROTEIN_FASTA_REGISTERED_BUT_NOT_LOCAL": "Download the registered protein FASTA, then rerun protein QC.",
        "PROTEIN_QC_NOT_AVAILABLE": "Run protein QC before using this resource as OrthoFinder input.",
    }
    return actions.get(blocker, "Manual review required.")


def _orthofinder_notes(scope: str, protein_row: dict[str, str] | None) -> str:
    if scope == "CDS_SUBSET_NOT_GENOME_WIDE":
        return "Excluded from OrthoFinder because it is a CDS-derived subset, not a genome-wide proteome."
    if scope == "NO_PROTEIN_FASTA":
        return "No genome-wide protein FASTA was found in the registered public sources."
    if protein_row:
        return protein_row.get("notes", "")
    return "NOT_ASSESSED"


def parse_gff_gene_coordinates(gff_path: Path, resource_id: str, accession: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with gff_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"{gff_path}:{line_number} does not have 9 GFF fields")
            seqid, _source, feature_type, start, end, _score, strand, _phase, attributes = fields
            if feature_type != "gene":
                continue
            parsed_attributes = _parse_gff_attributes(attributes)
            gene_id = parsed_attributes.get("ID") or parsed_attributes.get("gene_id") or "NOT_ASSESSED"
            gene_symbol = (
                parsed_attributes.get("Name")
                or parsed_attributes.get("gene")
                or parsed_attributes.get("gene_name")
                or "NOT_ASSESSED"
            )
            rows.append(
                {
                    "resource_id": resource_id,
                    "accession": accession,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "source_file": str(gff_path),
                    "parse_status": "PARSED_FROM_GFF",
                    "notes": "Coordinates parsed from local GFF. Coordinate parsing is not biological validation.",
                }
            )
    return rows


def _parse_gff_attributes(attributes: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in attributes.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        parsed[key.strip()] = value.strip().strip('"')
    return parsed


def build_gene_coordinate_rows(inventory_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    resource_ids = sorted({row["resource_id"] for row in inventory_rows})
    rows: list[dict[str, str]] = []
    for resource_id in resource_ids:
        gff_row = _first_by_role(inventory_rows, resource_id, "annotation_gff")
        if gff_row and gff_row["file_status"] == "LOCAL_AVAILABLE":
            rows.extend(
                parse_gff_gene_coordinates(Path(gff_row["local_path"]), resource_id, gff_row["accession"])
            )
            continue
        rows.append(
            {
                "resource_id": resource_id,
                "accession": gff_row["accession"] if gff_row else "NOT_ASSESSED",
                "gene_id": "NOT_ASSESSED",
                "gene_symbol": "NOT_ASSESSED",
                "seqid": "NOT_ASSESSED",
                "start": "NOT_ASSESSED",
                "end": "NOT_ASSESSED",
                "strand": "NOT_ASSESSED",
                "source_file": gff_row["local_path"] if gff_row else "NOT_AVAILABLE",
                "parse_status": "NO_LOCAL_GFF_OR_GTF",
                "notes": "No local GFF/GTF was available from searched public sources; do not infer gene absence.",
            }
        )
    return rows


def build_phase3_outputs(
    config_path: Path,
    inventory_path: Path,
    protein_qc_path: Path,
    source_checks_output: Path,
    gene_coordinates_output: Path,
    orthofinder_input_output: Path,
    copy_orthofinder_inputs: bool = False,
) -> dict[str, Path]:
    config = read_yaml(config_path)
    inventory_rows = read_tsv(inventory_path, REFERENCE_FILE_COLUMNS)
    protein_qc_rows = read_tsv(protein_qc_path, REFERENCE_PROTEIN_QC_COLUMNS)

    source_rows = configured_source_checks(config)
    coordinate_rows = build_gene_coordinate_rows(inventory_rows)
    orthofinder_rows = build_orthofinder_input_manifest(
        config,
        inventory_rows,
        protein_qc_rows,
        copy_inputs=copy_orthofinder_inputs,
    )

    write_tsv(source_checks_output, source_rows, PHASE3_SOURCE_CHECK_COLUMNS)
    write_tsv(gene_coordinates_output, coordinate_rows, REFERENCE_GENE_COORDINATE_COLUMNS)
    write_tsv(orthofinder_input_output, orthofinder_rows, ORTHOFINDER_INPUT_MANIFEST_COLUMNS)
    return {
        "source_checks": source_checks_output,
        "gene_coordinates": gene_coordinates_output,
        "orthofinder_input": orthofinder_input_output,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Phase 3 input discovery and OrthoFinder readiness outputs.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--protein-qc", type=Path, required=True)
    parser.add_argument("--source-checks-output", type=Path, required=True)
    parser.add_argument("--gene-coordinates-output", type=Path, required=True)
    parser.add_argument("--orthofinder-input-output", type=Path, required=True)
    parser.add_argument("--copy-orthofinder-inputs", action="store_true")
    args = parser.parse_args()

    configure_logging()
    build_phase3_outputs(
        args.config,
        args.inventory,
        args.protein_qc,
        args.source_checks_output,
        args.gene_coordinates_output,
        args.orthofinder_input_output,
        copy_orthofinder_inputs=args.copy_orthofinder_inputs,
    )


if __name__ == "__main__":
    main()
