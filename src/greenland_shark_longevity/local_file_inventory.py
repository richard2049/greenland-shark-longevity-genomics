"""Inventory manually downloaded public files."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from .fasta import parse_fasta
from .schemas import LOCAL_FILE_INVENTORY_COLUMNS
from .utils import configure_logging, write_tsv


def checksum(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest().upper()


def detect_file_role(path: Path) -> tuple[str, str, str]:
    name = path.name.lower()
    suffix = path.suffix.lower()
    if name.endswith(".gff") or name.endswith(".gff3"):
        return "annotation_gff", "gff3", "Parsed gene coordinates if registered as annotation_gff."
    if name.endswith(".gtf"):
        return "annotation_gtf", "gtf", "GTF coordinate parsing is not yet implemented."
    if name == "complete.proteins.faa" or suffix in {".faa", ".fa", ".fasta"} and "protein" in name:
        return "genome_annotation_protein_fasta", "protein_fasta", "Candidate genome-annotation protein set for QC and OrthoFinder staging."
    if name.endswith(".pep"):
        return "transcriptome_transdecoder_peptide_fasta", "protein_fasta", "Transcriptome-derived peptide predictions; keep separate from genome annotation proteome."
    if name.endswith(".fna") or name.endswith(".ffn"):
        return "coding_sequence_fasta", "nucleotide_fasta", "Coding-sequence FASTA; not a protein FASTA."
    if name.endswith(".zip"):
        return "source_archive", "zip", "Downloaded source archive; preserve for provenance."
    return "manual_review_required", suffix.lstrip(".") or "unknown", "File type requires manual review."


def sequence_count(path: Path, detected_format: str) -> str:
    if detected_format in {"protein_fasta", "nucleotide_fasta"}:
        records, _duplicates = parse_fasta(path)
        return str(len(records))
    return "NOT_APPLICABLE"


def inventory_local_files(root: Path, resource_id: str, source_id: str) -> list[dict[str, str]]:
    if not root.exists():
        raise FileNotFoundError(f"Manual file directory does not exist: {root}")
    rows: list[dict[str, str]] = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        role, detected_format, notes = detect_file_role(path)
        rows.append(
            {
                "resource_id": resource_id,
                "source_id": source_id,
                "file_name": path.name,
                "file_role_candidate": role,
                "detected_format": detected_format,
                "local_path": str(path),
                "byte_size": str(path.stat().st_size),
                "md5": checksum(path, "md5"),
                "sha256": checksum(path, "sha256"),
                "sequence_count": sequence_count(path, detected_format),
                "notes": notes,
            }
        )
    if not rows:
        raise ValueError(f"No files found under {root}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Inventory manually downloaded files with sizes and checksums.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--resource-id", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    configure_logging()
    rows = inventory_local_files(args.input_dir, args.resource_id, args.source_id)
    write_tsv(args.output, rows, LOCAL_FILE_INVENTORY_COLUMNS)


if __name__ == "__main__":
    main()
