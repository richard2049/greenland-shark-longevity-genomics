"""FASTA parsing and basic QC metrics."""

from __future__ import annotations

import argparse
import gzip
import logging
from pathlib import Path
from statistics import mean, median
from typing import TextIO

from .schemas import ASSEMBLY_QC_COLUMNS, PROTEIN_QC_COLUMNS
from .utils import configure_logging, write_tsv


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_fasta(path: Path) -> tuple[list[tuple[str, str]], list[str]]:
    records: list[tuple[str, str]] = []
    duplicate_ids: list[str] = []
    seen: set[str] = set()
    current_id: str | None = None
    current_parts: list[str] = []

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_parts).upper()))
                current_id = line[1:].split()[0]
                if not current_id:
                    raise ValueError(f"Empty FASTA identifier in {path} line {line_number}")
                if current_id in seen:
                    duplicate_ids.append(current_id)
                seen.add(current_id)
                current_parts = []
            else:
                if current_id is None:
                    raise ValueError(f"Sequence encountered before first header in {path} line {line_number}")
                current_parts.append(line)

    if current_id is not None:
        records.append((current_id, "".join(current_parts).upper()))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records, duplicate_ids


def n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    half = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= half:
            return length
    return 0


def assembly_qc(path: Path, resource_id: str) -> dict[str, str | int | float]:
    records, _duplicates = parse_fasta(path)
    sequences = [sequence for _record_id, sequence in records]
    lengths = [len(sequence) for sequence in sequences]
    total_length = sum(lengths)
    gc_count = sum(sequence.count("G") + sequence.count("C") for sequence in sequences)
    n_count = sum(sequence.count("N") for sequence in sequences)
    gc_percent = round((gc_count / total_length) * 100, 3) if total_length else 0.0
    return {
        "resource_id": resource_id,
        "sequence_count": len(records),
        "total_length_bp": total_length,
        "n50_bp": n50(lengths),
        "gc_percent": gc_percent,
        "n_count": n_count,
        "longest_sequence_bp": max(lengths),
        "busco_lineage": "NOT_ASSESSED",
        "busco_complete_percent": "NOT_ASSESSED",
        "busco_single_copy_percent": "NOT_ASSESSED",
        "busco_duplicated_percent": "NOT_ASSESSED",
        "busco_fragmented_percent": "NOT_ASSESSED",
        "busco_missing_percent": "NOT_ASSESSED",
        "notes": "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE" if "demo" in resource_id.lower() else "BUSCO not run by MVP QC script.",
    }


def protein_qc(
    path: Path,
    resource_id: str,
    ambiguous_residues: set[str] | None = None,
) -> dict[str, str | int | float]:
    ambiguous_residues = ambiguous_residues or {"B", "J", "O", "U", "X", "Z"}
    records, duplicate_ids = parse_fasta(path)
    sequences = [sequence for _record_id, sequence in records]
    lengths = [len(sequence.rstrip("*")) for sequence in sequences]
    ambiguous_count = sum(sum(1 for aa in sequence if aa in ambiguous_residues) for sequence in sequences)
    stop_count = sum(1 for sequence in sequences if "*" in sequence)
    return {
        "resource_id": resource_id,
        "protein_count": len(records),
        "total_aa": sum(lengths),
        "mean_length_aa": round(mean(lengths), 3),
        "median_length_aa": round(median(lengths), 3),
        "duplicate_ids": ",".join(sorted(set(duplicate_ids))) if duplicate_ids else "NONE",
        "ambiguous_residue_count": ambiguous_count,
        "sequences_with_stop_codon": stop_count,
        "notes": "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE" if "demo" in resource_id.lower() else "Protein QC only; no domain inference.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run basic FASTA QC.")
    parser.add_argument("--kind", choices=["assembly", "protein"], required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resource-id", required=True)
    args = parser.parse_args()

    configure_logging()
    logging.info("Running %s QC for %s", args.kind, args.input)
    if args.kind == "assembly":
        write_tsv(args.output, [assembly_qc(args.input, args.resource_id)], ASSEMBLY_QC_COLUMNS)
    else:
        write_tsv(args.output, [protein_qc(args.input, args.resource_id)], PROTEIN_QC_COLUMNS)


if __name__ == "__main__":
    main()
