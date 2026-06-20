"""Reference resource intake and Phase 1 resource-quality QC."""

from __future__ import annotations

import argparse
import hashlib
import logging
import urllib.request
from pathlib import Path

from .fasta import protein_qc
from .schemas import (
    REFERENCE_ANNOTATION_QC_COLUMNS,
    REFERENCE_ASSEMBLY_QC_COLUMNS,
    REFERENCE_FILE_COLUMNS,
    REFERENCE_PROTEIN_QC_COLUMNS,
)
from .utils import as_bool, configure_logging, read_tsv, read_yaml, write_tsv

NOT_AVAILABLE = {"", "TODO", "NOT_ASSESSED", "NOT_AVAILABLE"}


def _as_text(value: object) -> str:
    if value is None:
        return ""
    return str(value)


def _checksum(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest().upper()


def configured_reference_files(config: dict) -> list[dict[str, str]]:
    retrieval_date = _as_text(config.get("reference_intake", {}).get("retrieval_date", "NOT_ASSESSED"))
    rows = []
    for row in config.get("reference_files", []):
        normalized = {column: _as_text(row.get(column, "")) for column in REFERENCE_FILE_COLUMNS}
        normalized["selected_for_download"] = str(as_bool(row.get("selected_for_download")))
        normalized["download_method"] = normalized["download_method"] or "direct_url"
        normalized["retrieval_date"] = normalized["retrieval_date"] or retrieval_date
        rows.append(normalized)
    if not rows:
        raise ValueError("Config must include reference_files for reference intake")
    return rows


def download_selected(config_path: Path, inventory_output: Path) -> list[dict[str, str]]:
    rows = configured_reference_files(read_yaml(config_path))
    inventory_rows: list[dict[str, str]] = []
    for row in rows:
        url = row["url"]
        local_path = row["local_path"]
        selected = as_bool(row["selected_for_download"])
        download_method = row.get("download_method", "direct_url")

        if url in NOT_AVAILABLE:
            row["file_status"] = "NOT_AVAILABLE_FROM_REGISTERED_SOURCE"
            row["byte_size"] = "NOT_ASSESSED"
        elif selected:
            if local_path in NOT_AVAILABLE:
                raise ValueError(f"{row['resource_id']} {row['file_role']} selected for download without local_path")
            path = Path(local_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists():
                if download_method == "manual_browser":
                    row["file_status"] = "MISSING_LOCAL_MANUAL_DOWNLOAD_REQUIRED"
                    row["byte_size"] = "0"
                    inventory_rows.append(row)
                    continue
                logging.info("Downloading %s to %s", url, path)
                with urllib.request.urlopen(url, timeout=120) as response:
                    path.write_bytes(response.read())
            row["file_status"] = "LOCAL_AVAILABLE" if path.exists() else "MISSING_LOCAL"
            row["byte_size"] = str(path.stat().st_size) if path.exists() else "0"
            if path.exists():
                row["md5"] = row.get("md5") or _checksum(path, "md5")
                row["sha256"] = row.get("sha256") or _checksum(path, "sha256")
        else:
            row["file_status"] = "REMOTE_REGISTERED_NOT_DOWNLOADED"
            row["byte_size"] = "NOT_ASSESSED"
        inventory_rows.append(row)

    write_tsv(inventory_output, inventory_rows, REFERENCE_FILE_COLUMNS)
    return inventory_rows


def parse_assembly_stats(path: Path) -> tuple[dict[str, str], dict[str, str], str]:
    metadata: dict[str, str] = {}
    stats_by_scope: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                stripped = line.lstrip("#").strip()
                if ":" in stripped and not stripped.startswith("Statistic"):
                    key, value = stripped.split(":", 1)
                    metadata[key.strip().lower()] = value.strip()
                continue
            fields = line.split("\t")
            if len(fields) != 6:
                continue
            unit_name, molecule_name, molecule_type, sequence_type, statistic, value = fields
            stats_by_scope.setdefault((unit_name, molecule_name, molecule_type, sequence_type), {})[statistic] = value

    preferred_key = ("Primary Assembly", "all", "all", "all")
    fallback_key = ("all", "all", "all", "all")
    if preferred_key in stats_by_scope:
        merged = dict(stats_by_scope.get(fallback_key, {}))
        merged.update(stats_by_scope[preferred_key])
        return metadata, merged, "Primary Assembly"
    if fallback_key in stats_by_scope:
        return metadata, stats_by_scope[fallback_key], "all"
    raise ValueError(f"Could not find primary/all assembly statistics in {path}")


def parse_feature_count(path: Path) -> tuple[int, int]:
    with path.open("r", encoding="utf-8") as handle:
        header: list[str] | None = None
        rows = []
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if header is None:
                header = line.lstrip("# ").split("\t")
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            rows.append(dict(zip(header, fields)))
    if header is None:
        raise ValueError(f"{path} has no feature-count header")
    required = {"Feature", "Class", "Assembly-unit name", "Unique Ids"}
    missing = required - set(header)
    if missing:
        raise ValueError(f"{path} is missing feature-count columns: {', '.join(sorted(missing))}")
    scoped_rows = [row for row in rows if row["Assembly-unit name"] == "all"] or rows
    protein_coding = 0
    total_gene = 0
    for row in scoped_rows:
        if row["Feature"] != "gene":
            continue
        value = int(row["Unique Ids"])
        total_gene += value
        if row["Class"] == "protein_coding":
            protein_coding += value
    return protein_coding, total_gene


def parse_gff_gene_counts(path: Path) -> tuple[int, int]:
    protein_coding = 0
    total_gene = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number} does not have 9 GFF fields")
            if fields[2] != "gene":
                continue
            total_gene += 1
            attributes = fields[8]
            if "gene_biotype=protein_coding" in attributes or "gbkey=Gene" in attributes:
                protein_coding += 1
    return protein_coding, total_gene


def _first_by_role(rows: list[dict[str, str]], resource_id: str, role: str) -> dict[str, str] | None:
    for row in rows:
        if row["resource_id"] == resource_id and row["file_role"] == role:
            return row
    return None


def build_reference_qc(config_path: Path, inventory_path: Path, output_dir: Path) -> dict[str, Path]:
    _config = read_yaml(config_path)
    inventory = read_tsv(inventory_path, REFERENCE_FILE_COLUMNS)
    resource_ids = sorted({row["resource_id"] for row in inventory})
    assembly_rows: list[dict[str, str]] = []
    annotation_rows: list[dict[str, str]] = []
    protein_rows: list[dict[str, str]] = []

    for resource_id in resource_ids:
        stats_row = _first_by_role(inventory, resource_id, "assembly_stats")
        feature_row = _first_by_role(inventory, resource_id, "feature_count")
        protein_row = _first_by_role(inventory, resource_id, "protein_fasta")
        gff_row = _first_by_role(inventory, resource_id, "annotation_gff")

        accession = stats_row["accession"] if stats_row else (protein_row["accession"] if protein_row else resource_id)
        if stats_row and stats_row["file_status"] == "LOCAL_AVAILABLE":
            metadata, stats, scope = parse_assembly_stats(Path(stats_row["local_path"]))
            assembly_rows.append(
                {
                    "resource_id": resource_id,
                    "accession": accession,
                    "assembly_name": metadata.get("assembly name", "NOT_ASSESSED"),
                    "bioproject": metadata.get("bioproject", "NOT_ASSESSED"),
                    "biosample": metadata.get("biosample", "NOT_ASSESSED"),
                    "assembly_level": metadata.get("assembly level", "NOT_ASSESSED"),
                    "release_date": metadata.get("date", "NOT_ASSESSED"),
                    "source_file": stats_row["local_path"],
                    "sequence_scope": scope,
                    "total_length_bp": stats.get("total-length", "NOT_ASSESSED"),
                    "ungapped_length_bp": stats.get("ungapped-length", "NOT_ASSESSED"),
                    "scaffold_count": stats.get("scaffold-count", "NOT_ASSESSED"),
                    "scaffold_n50_bp": stats.get("scaffold-N50", "NOT_ASSESSED"),
                    "scaffold_l50": stats.get("scaffold-L50", "NOT_ASSESSED"),
                    "contig_count": stats.get("contig-count", "NOT_ASSESSED"),
                    "contig_n50_bp": stats.get("contig-N50", "NOT_ASSESSED"),
                    "gc_percent": stats.get("gc-perc", "NOT_ASSESSED"),
                    "molecule_count": stats.get("molecule-count", "NOT_ASSESSED"),
                    "top_level_count": stats.get("top-level-count", "NOT_ASSESSED"),
                    "busco_lineage": "NOT_ASSESSED",
                    "busco_complete_percent": "NOT_ASSESSED",
                    "notes": "Resource-quality QC from NCBI assembly_stats.txt; not biological evidence.",
                }
            )

        if gff_row and gff_row["file_status"] == "LOCAL_AVAILABLE":
            protein_coding, total_gene = parse_gff_gene_counts(Path(gff_row["local_path"]))
            annotation_rows.append(
                {
                    "resource_id": resource_id,
                    "accession": gff_row["accession"],
                    "source_file": gff_row["local_path"],
                    "protein_coding_gene_count": protein_coding,
                    "total_gene_count": total_gene,
                    "separate_gff_url": gff_row["url"],
                    "separate_gff_status": gff_row["file_status"],
                    "notes": "Gene counts parsed from local Figshare GFF. This is annotation resource QC, not biological evidence for candidate mechanisms.",
                }
            )
        elif feature_row and feature_row["file_status"] == "LOCAL_AVAILABLE":
            protein_coding, total_gene = parse_feature_count(Path(feature_row["local_path"]))
            annotation_rows.append(
                {
                    "resource_id": resource_id,
                    "accession": feature_row["accession"],
                    "source_file": feature_row["local_path"],
                    "protein_coding_gene_count": protein_coding,
                    "total_gene_count": total_gene,
                    "separate_gff_url": gff_row["url"] if gff_row else "NOT_ASSESSED",
                    "separate_gff_status": gff_row["file_status"] if gff_row else "NOT_ASSESSED",
                    "notes": "Feature count reports annotation availability only; zero protein-coding genes means no protein-coding annotation in this NCBI assembly package, not gene absence from the genome.",
                }
            )

        if protein_row is None or protein_row["url"] in NOT_AVAILABLE:
            protein_rows.append(
                {
                    "resource_id": resource_id,
                    "accession": accession,
                    "protein_fasta_url": "NOT_AVAILABLE",
                    "local_path": "NOT_AVAILABLE",
                    "status": "NO_SEPARATE_PROTEIN_FASTA_IN_REGISTERED_NCBI_ASSEMBLY_DIRECTORY",
                    "protein_count": "NOT_ASSESSED",
                    "total_aa": "NOT_ASSESSED",
                    "mean_length_aa": "NOT_ASSESSED",
                    "median_length_aa": "NOT_ASSESSED",
                    "duplicate_ids": "NOT_ASSESSED",
                    "ambiguous_residue_count": "NOT_ASSESSED",
                    "sequences_with_stop_codon": "NOT_ASSESSED",
                    "notes": "Protein sequence QC was not run because no separate protein FASTA is available for this registered assembly source.",
                }
            )
        elif protein_row["file_status"] == "LOCAL_AVAILABLE":
            qc = protein_qc(Path(protein_row["local_path"]), accession)
            protein_rows.append(
                {
                    "resource_id": resource_id,
                    "accession": accession,
                    "protein_fasta_url": protein_row["url"],
                    "local_path": protein_row["local_path"],
                    "status": "LOCAL_AVAILABLE",
                    "protein_count": str(qc["protein_count"]),
                    "total_aa": str(qc["total_aa"]),
                    "mean_length_aa": str(qc["mean_length_aa"]),
                    "median_length_aa": str(qc["median_length_aa"]),
                    "duplicate_ids": str(qc["duplicate_ids"]),
                    "ambiguous_residue_count": str(qc["ambiguous_residue_count"]),
                    "sequences_with_stop_codon": str(qc["sequences_with_stop_codon"]),
                    "notes": f"{protein_row['notes']} Basic FASTA QC only; no domain inference.",
                }
            )
        else:
            protein_rows.append(
                {
                    "resource_id": resource_id,
                    "accession": accession,
                    "protein_fasta_url": protein_row["url"],
                    "local_path": protein_row["local_path"],
                    "status": protein_row["file_status"],
                    "protein_count": "NOT_ASSESSED",
                    "total_aa": "NOT_ASSESSED",
                    "mean_length_aa": "NOT_ASSESSED",
                    "median_length_aa": "NOT_ASSESSED",
                    "duplicate_ids": "NOT_ASSESSED",
                    "ambiguous_residue_count": "NOT_ASSESSED",
                    "sequences_with_stop_codon": "NOT_ASSESSED",
                    "notes": "Protein FASTA is registered but has not been downloaded.",
                }
            )

    paths = {
        "assembly": output_dir / "reference_assembly_qc.tsv",
        "annotation": output_dir / "reference_annotation_qc.tsv",
        "protein": output_dir / "reference_protein_qc.tsv",
    }
    write_tsv(paths["assembly"], assembly_rows, REFERENCE_ASSEMBLY_QC_COLUMNS)
    write_tsv(paths["annotation"], annotation_rows, REFERENCE_ANNOTATION_QC_COLUMNS)
    write_tsv(paths["protein"], protein_rows, REFERENCE_PROTEIN_QC_COLUMNS)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description="Download selected reference files and run Phase 1 resource QC.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    download_parser = subparsers.add_parser("download")
    download_parser.add_argument("--config", type=Path, required=True)
    download_parser.add_argument("--inventory-output", type=Path, required=True)

    qc_parser = subparsers.add_parser("qc")
    qc_parser.add_argument("--config", type=Path, required=True)
    qc_parser.add_argument("--inventory", type=Path, required=True)
    qc_parser.add_argument("--output-dir", type=Path, required=True)

    args = parser.parse_args()
    configure_logging()
    if args.command == "download":
        download_selected(args.config, args.inventory_output)
    elif args.command == "qc":
        build_reference_qc(args.config, args.inventory, args.output_dir)


if __name__ == "__main__":
    main()
