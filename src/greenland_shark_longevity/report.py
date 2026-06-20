"""Generate a preliminary human-readable report."""

from __future__ import annotations

import argparse
from pathlib import Path

from .utils import ensure_parent, read_tsv


def _count_by(rows: list[dict[str, str]], column: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        key = row.get(column, "NOT_ASSESSED")
        counts[key] = counts.get(key, 0) + 1
    return counts


def generate_report(
    data_manifest: Path,
    species_manifest: Path,
    assembly_qc: Path,
    protein_qc: Path,
    evidence: Path,
    output: Path,
) -> None:
    data_rows = read_tsv(data_manifest)
    species_rows = read_tsv(species_manifest)
    assembly_rows = read_tsv(assembly_qc)
    protein_rows = read_tsv(protein_qc)
    evidence_rows = read_tsv(evidence)
    tier_counts = _count_by(evidence_rows, "evidence_tier")

    lines = [
        "# Preliminary MVP Report",
        "",
        "This report is generated from the demo workflow. Demo outputs are not biological evidence.",
        "",
        "## Registered Resources",
        "",
        f"- Resources registered: {len(data_rows)}",
        f"- Species entries: {len(species_rows)}",
        "",
        "## Demo QC",
        "",
    ]
    for row in assembly_rows:
        lines.append(
            f"- Assembly `{row['resource_id']}`: {row['sequence_count']} sequences, "
            f"{row['total_length_bp']} bp total, N50 {row['n50_bp']} bp, GC {row['gc_percent']}%."
        )
    for row in protein_rows:
        lines.append(
            f"- Protein set `{row['resource_id']}`: {row['protein_count']} proteins, "
            f"median length {row['median_length_aa']} aa, duplicate IDs {row['duplicate_ids']}."
        )
    lines.extend(["", "## Evidence Tier Summary", ""])
    for tier, count in sorted(tier_counts.items()):
        lines.append(f"- {tier}: {count}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The MVP currently answers workflow-readiness questions only: schemas validate, toy fixtures flow through QC and evidence scoring, and public resources are tracked before download.",
            "",
            "No robust biological mechanism is claimed from the demo run.",
            "",
            "## Next Validation Needed",
            "",
            "- Select exact public genome/protein/annotation files and record retrieval metadata.",
            "- Add real OrthoFinder or equivalent orthology outputs.",
            "- Validate candidate duplications with protein domains, coordinates, isoform filtering, and cross-resource checks.",
            "- Process retinal RNA-seq only after sample metadata and design are reviewed.",
        ]
    )
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate preliminary MVP report.")
    parser.add_argument("--data-manifest", type=Path, required=True)
    parser.add_argument("--species-manifest", type=Path, required=True)
    parser.add_argument("--assembly-qc", type=Path, required=True)
    parser.add_argument("--protein-qc", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    generate_report(
        args.data_manifest,
        args.species_manifest,
        args.assembly_qc,
        args.protein_qc,
        args.evidence,
        args.output,
    )


if __name__ == "__main__":
    main()

