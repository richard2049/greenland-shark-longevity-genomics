"""Candidate-specific validation readiness after OrthoFinder mapping."""

from __future__ import annotations

import argparse
from pathlib import Path

from .candidate_orthofinder import TARGET_RESOURCE_ID, candidate_aliases
from .candidate_panels import iter_candidates
from .schemas import (
    CANDIDATE_ISOFORM_AUDIT_COLUMNS,
    COPY_NUMBER_COLUMNS,
    DOMAIN_CHECK_TARGET_COLUMNS,
    DUPLICATION_AUDIT_COLUMNS,
    RESCUE_TARGET_COLUMNS,
)
from .utils import configure_logging, read_tsv, read_yaml, write_tsv


DEFAULT_HIGH_PRIORITY_RESCUE_GENES = ["H1F0", "FTH1B", "TP53", "RAD51"]


def split_list(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip() and part.strip() != "NOT_ASSESSED"]


def parse_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current_id: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        sequences[current_id] = "".join(chunks)
    return sequences


def protein_variants(protein_id: str) -> list[str]:
    return [
        protein_id,
        protein_id.replace("_", ":").replace("WGS:ZZZZ", "WGS:ZZZZ"),
        protein_id.replace("|WGS_ZZZZ|", "|WGS:ZZZZ|"),
        protein_id.replace("|", "_").replace(":", "_"),
    ]


def find_sequence(protein_id: str, sequences: dict[str, str]) -> tuple[str, str]:
    for variant in protein_variants(protein_id):
        if variant in sequences:
            return variant, sequences[variant]
    return protein_id, ""


def representative_for_gene(protein_ids: list[str], sequences: dict[str, str]) -> tuple[str, str, int, str, str]:
    observed: list[tuple[str, str, int]] = []
    for protein_id in protein_ids:
        fasta_id, sequence = find_sequence(protein_id, sequences)
        observed.append((protein_id, fasta_id, len(sequence)))
    available = [item for item in observed if item[2] > 0]
    if not available:
        return "NOT_ASSESSED", "NOT_ASSESSED", 0, "sequence_missing", "No candidate protein sequence was found."
    representative = max(available, key=lambda item: (item[2], item[0]))
    return representative[0], representative[1], representative[2], "longest_protein_isoform", "Representative selected as longest available protein isoform for the locus."


def candidate_lookup(panel_path: Path) -> dict[str, dict[str, object]]:
    return {str(candidate["gene_symbol"]).upper(): candidate for candidate in iter_candidates(read_yaml(panel_path))}


def build_isoform_and_domain_rows(
    candidate_panels: Path,
    candidate_copy_number: Path,
    duplication_audit: Path,
    protein_fasta: Path,
    representative_fasta: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, str]]:
    copy_rows = read_tsv(candidate_copy_number, COPY_NUMBER_COLUMNS)
    audit_rows = read_tsv(duplication_audit, DUPLICATION_AUDIT_COLUMNS)
    audit_by_gene = {row["gene_symbol"]: row for row in audit_rows}
    sequences = parse_fasta(protein_fasta)

    isoform_rows: list[dict[str, str]] = []
    domain_rows: list[dict[str, str]] = []
    representative_sequences: dict[str, str] = {}

    for copy_row in copy_rows:
        gene_symbol = copy_row["gene_symbol"]
        protein_ids = split_list(copy_row["protein_ids"])
        gene_ids = split_list(copy_row["gene_ids"])
        if not gene_ids or not protein_ids:
            continue
        audit_row = audit_by_gene.get(gene_symbol, {})
        proteins_by_gene: dict[str, list[str]] = {gene_id: [] for gene_id in gene_ids}
        if len(gene_ids) == 1:
            proteins_by_gene[gene_ids[0]] = protein_ids
        else:
            for protein_id in protein_ids:
                matched_gene = next((gene_id for gene_id in gene_ids if gene_id.replace("gene-", "") in protein_id), "")
                if matched_gene:
                    proteins_by_gene[matched_gene].append(protein_id)

        for gene_id, locus_proteins in proteins_by_gene.items():
            representative_id, fasta_id, rep_length, method, note = representative_for_gene(locus_proteins, sequences)
            lengths: list[str] = []
            for protein_id in locus_proteins:
                _fasta_id, sequence = find_sequence(protein_id, sequences)
                lengths.append(str(len(sequence)) if sequence else "MISSING")
            isoform_count = len(locus_proteins)
            status = "single_protein" if isoform_count == 1 else "representative_selected"
            if representative_id == "NOT_ASSESSED":
                status = "sequence_missing"
            representative_fasta_id = f"{gene_symbol}|{gene_id}|{representative_id}"
            if fasta_id != "NOT_ASSESSED" and rep_length > 0:
                representative_sequences[representative_fasta_id] = sequences[fasta_id]

            isoform_rows.append(
                {
                    "mechanism": copy_row["mechanism"],
                    "gene_symbol": gene_symbol,
                    "resource_id": copy_row["resource_id"],
                    "gene_id": gene_id,
                    "orthogroup_id": copy_row["orthogroup_id"],
                    "coordinate_summary": audit_row.get("coordinate_summary", "NOT_ASSESSED"),
                    "protein_ids": ",".join(locus_proteins) if locus_proteins else "NOT_ASSESSED",
                    "protein_lengths_aa": ",".join(lengths) if lengths else "NOT_ASSESSED",
                    "isoform_count": str(isoform_count),
                    "representative_protein_id": representative_id,
                    "representative_length_aa": str(rep_length) if rep_length else "NOT_ASSESSED",
                    "representative_fasta_id": representative_fasta_id if rep_length else "NOT_ASSESSED",
                    "selection_method": method,
                    "isoform_filter_status": status,
                    "notes": note,
                }
            )

            domain_rows.append(
                {
                    "mechanism": copy_row["mechanism"],
                    "gene_symbol": gene_symbol,
                    "resource_id": copy_row["resource_id"],
                    "gene_id": gene_id,
                    "orthogroup_id": copy_row["orthogroup_id"],
                    "representative_protein_id": representative_id,
                    "representative_fasta_id": representative_fasta_id if rep_length else "NOT_ASSESSED",
                    "representative_length_aa": str(rep_length) if rep_length else "NOT_ASSESSED",
                    "candidate_fasta": str(representative_fasta),
                    "domain_check_status": "READY_FOR_DOMAIN_SCAN" if rep_length else "SEQUENCE_UNAVAILABLE",
                    "recommended_method": "Run InterProScan or HMMER/Pfam on representative candidate proteins; do not infer domain integrity from product names.",
                    "required_validation": "REQUIRES_VALIDATION: record database version, domain accession(s), coverage, e-value, and whether expected domains are complete.",
                    "notes": "Domain integrity is not assessed by this repository step.",
                }
            )

    return isoform_rows, domain_rows, representative_sequences


def write_representative_fasta(path: Path, sequences: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for fasta_id, sequence in sorted(sequences.items()):
            handle.write(f">{fasta_id}\n")
            for idx in range(0, len(sequence), 80):
                handle.write(sequence[idx : idx + 80] + "\n")


def build_rescue_rows(
    candidate_panels: Path,
    candidate_copy_number: Path,
    high_priority_genes: list[str] | None = None,
) -> list[dict[str, str]]:
    high_priority_genes = [gene.upper() for gene in (high_priority_genes or DEFAULT_HIGH_PRIORITY_RESCUE_GENES)]
    candidates = candidate_lookup(candidate_panels)
    copy_rows = {row["gene_symbol"]: row for row in read_tsv(candidate_copy_number, COPY_NUMBER_COLUMNS)}
    rows: list[dict[str, str]] = []
    for priority, gene_symbol in enumerate(high_priority_genes, start=1):
        candidate = candidates.get(gene_symbol)
        copy_row = copy_rows.get(gene_symbol, {})
        aliases = sorted(candidate_aliases(candidate)) if candidate else [gene_symbol]
        mapping_status = copy_row.get("mapping_status", "NOT_ASSESSED")
        if mapping_status == "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH":
            method = "Reciprocal protein similarity search against staged proteomes, then protein-to-genome alignment if genome FASTA is available."
            blockers = "Validated query/reference protein set and local genome FASTA are required for rescue; do not infer absence from current annotation."
        elif copy_row.get("protein_ids", "NOT_ASSESSED") == "NOT_ASSESSED":
            method = "Inspect GFF/CDS/protein ID consistency and rerun OrthoFinder mapping if needed."
            blockers = "Exact symbol exists but no mapped OrthoFinder protein was recorded."
        else:
            method = "Run representative-protein domain scan and inspect orthogroup membership before interpretation."
            blockers = "Domain database/tool output is not available yet."
        rows.append(
            {
                "priority": str(priority),
                "mechanism": str(candidate.get("mechanism", "NOT_ASSESSED")) if candidate else "NOT_ASSESSED",
                "gene_symbol": gene_symbol,
                "mapping_status": mapping_status,
                "current_copy_count": copy_row.get("copy_count", "NOT_ASSESSED"),
                "current_orthogroup_id": copy_row.get("orthogroup_id", "NOT_ASSESSED"),
                "aliases_considered": ",".join(aliases),
                "recommended_next_method": method,
                "required_inputs": "representative protein FASTA; domain database; genome FASTA for rescue; comparator proteomes for reciprocal checks",
                "blockers": blockers,
                "claim_language_guardrail": "Do not report absent, inactivated, activated, adaptive, or causal language from this target table.",
                "notes": "High-priority candidate validation target generated from current exact-symbol OrthoFinder mapping.",
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Build candidate-specific isoform, domain, and rescue validation targets.")
    parser.add_argument("--candidate-panels", type=Path, required=True)
    parser.add_argument("--candidate-copy-number", type=Path, required=True)
    parser.add_argument("--duplication-audit", type=Path, required=True)
    parser.add_argument("--protein-fasta", type=Path, required=True)
    parser.add_argument("--representative-fasta-output", type=Path, required=True)
    parser.add_argument("--isoform-audit-output", type=Path, required=True)
    parser.add_argument("--domain-targets-output", type=Path, required=True)
    parser.add_argument("--rescue-targets-output", type=Path, required=True)
    parser.add_argument("--high-priority-genes", nargs="*", default=DEFAULT_HIGH_PRIORITY_RESCUE_GENES)
    args = parser.parse_args()

    configure_logging()
    isoform_rows, domain_rows, representative_sequences = build_isoform_and_domain_rows(
        args.candidate_panels,
        args.candidate_copy_number,
        args.duplication_audit,
        args.protein_fasta,
        args.representative_fasta_output,
    )
    rescue_rows = build_rescue_rows(args.candidate_panels, args.candidate_copy_number, args.high_priority_genes)
    write_representative_fasta(args.representative_fasta_output, representative_sequences)
    write_tsv(args.isoform_audit_output, isoform_rows, CANDIDATE_ISOFORM_AUDIT_COLUMNS)
    write_tsv(args.domain_targets_output, domain_rows, DOMAIN_CHECK_TARGET_COLUMNS)
    write_tsv(args.rescue_targets_output, rescue_rows, RESCUE_TARGET_COLUMNS)


if __name__ == "__main__":
    main()
