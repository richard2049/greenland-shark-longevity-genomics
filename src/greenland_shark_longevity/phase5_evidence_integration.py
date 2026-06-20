"""Integrate Phase 5 repeat-context artifact risk into evidence outputs."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .evidence import validate_guardrail_language
from .schemas import EVIDENCE_COLUMNS, PHASE5_GENE_REPEAT_CONTEXT_COLUMNS, PHASE5C_GENE_QC_COLUMNS
from .utils import ensure_parent, read_tsv, write_tsv

LOGGER = logging.getLogger(__name__)

PHASE5_GENES = {"FTH1B", "H1F0", "RAD51", "TP53"}
NOT_ASSESSED = "NOT_ASSESSED"


def parse_int(value: str | None) -> int:
    if value in {None, "", NOT_ASSESSED}:
        return 0
    return int(float(str(value)))


def parse_float(value: str | None) -> float:
    if value in {None, "", NOT_ASSESSED}:
        return 0.0
    return float(str(value))


def append_unique_semicolon(value: str, additions: list[str]) -> str:
    parts: list[str] = []
    for item in [value, *additions]:
        if item in {"", NOT_ASSESSED}:
            continue
        for part in str(item).split(";"):
            cleaned = part.strip()
            if cleaned and cleaned not in parts:
                parts.append(cleaned)
    return ";".join(parts) if parts else NOT_ASSESSED


def repeat_artifact_label(summary_row: dict[str, str], base_risk: str) -> str:
    overlap_loci = parse_int(summary_row.get("loci_with_repeat_overlap"))
    max_fraction = parse_float(summary_row.get("max_locus_repeat_overlap_fraction"))
    status = summary_row.get("repeat_annotation_status", NOT_ASSESSED)

    if status != "REPEAT_ANNOTATION_AVAILABLE":
        return append_unique_semicolon(base_risk, ["repeat_context_not_assessed"])
    if overlap_loci > 0 and max_fraction >= 0.5:
        return append_unique_semicolon(base_risk, ["high_repeat_context_artifact_risk"])
    if overlap_loci > 0:
        return append_unique_semicolon(base_risk, ["moderate_repeat_context_artifact_risk"])
    return append_unique_semicolon(base_risk, ["repeat_context_window_only_no_direct_overlap"])


def repeat_context_sentence(summary_row: dict[str, str]) -> str:
    overlap_loci = parse_int(summary_row.get("loci_with_repeat_overlap"))
    locus_count = parse_int(summary_row.get("locus_count"))
    total_overlap = parse_int(summary_row.get("total_repeat_overlap_bp"))
    max_fraction = parse_float(summary_row.get("max_locus_repeat_overlap_fraction"))
    if summary_row.get("repeat_annotation_status") != "REPEAT_ANNOTATION_AVAILABLE":
        return (
            "Phase 5 repeat context is not assessable from the current local files and remains an artifact-risk gap."
        )
    if overlap_loci > 0:
        return (
            f"Phase 5 records repeat overlap in {overlap_loci}/{locus_count} reviewed loci "
            f"({total_overlap} bp total overlap; maximum locus fraction {max_fraction:.6g}) as artifact/context evidence only."
        )
    return (
        f"Phase 5 records repeat features in the local candidate window but no direct repeat overlap across {locus_count} reviewed loci."
    )


def build_phase5_evidence_row(
    preferred_row: dict[str, str],
    summary_row: dict[str, str],
    phase5c_qc_row: dict[str, str] | None = None,
) -> dict[str, str]:
    gene_symbol = summary_row["gene_symbol"]
    repeat_sentence = repeat_context_sentence(summary_row)
    supporting_files = summary_row.get("supporting_files", NOT_ASSESSED)
    existing_interpretation = preferred_row.get("biological_interpretation", NOT_ASSESSED)
    existing_validation = preferred_row.get("required_validation", NOT_ASSESSED)
    existing_guardrail = preferred_row.get("claim_language_guardrail", NOT_ASSESSED)

    phase5c_sentence = ""
    phase5c_supporting_files = NOT_ASSESSED
    phase5c_required_validation = NOT_ASSESSED
    phase5c_guardrail = NOT_ASSESSED
    phase5c_artifact = NOT_ASSESSED
    if phase5c_qc_row is not None:
        phase5c_sentence = (
            f" Phase 5c QC status is {phase5c_qc_row['phase5c_qc_status']}; "
            f"RepeatMasker .out direct-overlap loci are {phase5c_qc_row['repeatmasker_out_loci_with_direct_overlap']}/"
            f"{phase5c_qc_row['locus_count']} with maximum window density {phase5c_qc_row['max_window_repeat_density']}."
        )
        phase5c_supporting_files = phase5c_qc_row.get("supporting_files", NOT_ASSESSED)
        phase5c_required_validation = phase5c_qc_row.get("required_validation", NOT_ASSESSED)
        phase5c_guardrail = phase5c_qc_row.get("claim_language_guardrail", NOT_ASSESSED)
        phase5c_artifact = phase5c_qc_row.get("artifact_risk_modifier", NOT_ASSESSED)

    biological_interpretation = (
        f"{existing_interpretation} {repeat_sentence} "
        f"{phase5c_sentence} "
        "This repeat context does not support repeat-mediated mechanism, validated duplication, or longevity interpretation."
    )
    required_validation = append_unique_semicolon(
        existing_validation,
        [
            summary_row.get("required_validation", NOT_ASSESSED),
            phase5c_required_validation,
            "Manual review of repeat annotations, locus coordinates, paralog identity, and cross-resource support before Phase 8 scoring.",
        ],
    )
    guardrail = append_unique_semicolon(
        existing_guardrail,
        [
            summary_row.get("claim_language_guardrail", NOT_ASSESSED),
            phase5c_guardrail,
            "Use Phase 5 repeat context only as artifact risk in reports and evidence scoring.",
        ],
    )
    resources = append_unique_semicolon(
        preferred_row.get("resources_supporting", NOT_ASSESSED),
        [
            supporting_files,
            phase5c_supporting_files,
            "RepeatModeler_2.0.8",
            "RepeatMasker_4.2.3",
            "SMIC_TOKYO_GENOME_2025",
        ],
    )

    row = {
        "mechanism": preferred_row.get("mechanism") or summary_row["mechanism"],
        "gene_or_pathway": gene_symbol,
        "evidence_type": append_unique_semicolon(
            preferred_row.get("evidence_type", NOT_ASSESSED),
            ["phase5_repeat_context_artifact_risk"],
        ),
        "evidence_tier": preferred_row.get("evidence_tier", "Artifact/uncertain"),
        "resources_supporting": resources,
        "artifact_risk": append_unique_semicolon(
            repeat_artifact_label(summary_row, preferred_row.get("artifact_risk", NOT_ASSESSED)),
            [phase5c_artifact],
        ),
        "biological_interpretation": biological_interpretation,
        "relevance_to_aging_longevity": preferred_row.get("relevance_to_aging_longevity", NOT_ASSESSED),
        "translational_relevance": preferred_row.get("translational_relevance", NOT_ASSESSED),
        "required_validation": required_validation,
        "claim_language_guardrail": guardrail,
    }
    validate_guardrail_language([row])
    return row


def integrate_phase5_repeat_context(
    base_evidence_path: Path,
    phase4e_evidence_path: Path,
    phase5_gene_summary_path: Path,
    phase5c_gene_qc_path: Path | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    base_rows = read_tsv(base_evidence_path, EVIDENCE_COLUMNS)
    phase4e_rows = read_tsv(phase4e_evidence_path, EVIDENCE_COLUMNS)
    summary_rows = read_tsv(phase5_gene_summary_path, PHASE5_GENE_REPEAT_CONTEXT_COLUMNS)
    phase5c_by_gene: dict[str, dict[str, str]] = {}
    if phase5c_gene_qc_path is not None and phase5c_gene_qc_path.exists():
        phase5c_by_gene = {
            row["gene_symbol"]: row
            for row in read_tsv(phase5c_gene_qc_path, PHASE5C_GENE_QC_COLUMNS)
        }
    if not summary_rows:
        raise ValueError(f"{phase5_gene_summary_path} has no Phase 5 gene summary rows")

    base_by_gene = {row["gene_or_pathway"]: row for row in base_rows}
    phase4e_by_gene = {row["gene_or_pathway"]: row for row in phase4e_rows}
    phase5_evidence_rows: list[dict[str, str]] = []
    for summary_row in sorted(summary_rows, key=lambda row: row["gene_symbol"]):
        gene_symbol = summary_row["gene_symbol"]
        preferred_row = phase4e_by_gene.get(gene_symbol) or base_by_gene.get(gene_symbol)
        if preferred_row is None:
            preferred_row = {
                "mechanism": summary_row["mechanism"],
                "gene_or_pathway": gene_symbol,
                "evidence_type": "phase5_repeat_context_artifact_risk",
                "evidence_tier": "Artifact/uncertain",
                "resources_supporting": NOT_ASSESSED,
                "artifact_risk": NOT_ASSESSED,
                "biological_interpretation": "Candidate repeat context requires integration with upstream evidence before interpretation.",
                "relevance_to_aging_longevity": NOT_ASSESSED,
                "translational_relevance": NOT_ASSESSED,
                "required_validation": "REQUIRES_VALIDATION",
                "claim_language_guardrail": "Do not interpret repeat context as biological mechanism.",
            }
        phase5_evidence_rows.append(build_phase5_evidence_row(preferred_row, summary_row, phase5c_by_gene.get(gene_symbol)))

    phase5_genes = {row["gene_or_pathway"] for row in phase5_evidence_rows}
    retained_rows = [row for row in base_rows if row.get("gene_or_pathway") not in phase5_genes]
    integrated_rows = retained_rows + phase5_evidence_rows
    validate_guardrail_language(integrated_rows)
    return phase5_evidence_rows, integrated_rows


def write_phase5_report(
    phase5_evidence_rows: list[dict[str, str]],
    phase5_gene_summary_rows: list[dict[str, str]],
    phase5c_gene_qc_rows: list[dict[str, str]] | None,
    output: Path,
) -> None:
    summary_by_gene = {row["gene_symbol"]: row for row in phase5_gene_summary_rows}
    phase5c_by_gene = {row["gene_symbol"]: row for row in phase5c_gene_qc_rows or []}
    lines = [
        "# Phase 5 Repeat Context Artifact-Risk Report",
        "",
        "Phase 5 records local repeat context around hardened candidate loci. These rows are artifact/context evidence only and do not support repeat-mediated mechanism, validated duplication, adaptation, pathway activity, or longevity interpretation.",
        "",
        "| Gene | Evidence tier | Repeat overlap | Phase 5c QC | Artifact risk | Interpretation guardrail |",
        "|---|---|---:|---|---|---|",
    ]
    for row in sorted(phase5_evidence_rows, key=lambda item: item["gene_or_pathway"]):
        gene = row["gene_or_pathway"]
        summary = summary_by_gene.get(gene, {})
        overlap = (
            f"{summary.get('loci_with_repeat_overlap', NOT_ASSESSED)}/"
            f"{summary.get('locus_count', NOT_ASSESSED)} loci; "
            f"max fraction {summary.get('max_locus_repeat_overlap_fraction', NOT_ASSESSED)}"
        )
        phase5c_status = phase5c_by_gene.get(gene, {}).get("phase5c_qc_status", "NOT_ASSESSED")
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{gene}`",
                    row["evidence_tier"],
                    overlap,
                    phase5c_status,
                    row["artifact_risk"],
                    "Repeat context is artifact risk only.",
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "Supporting tables:",
            "",
            "- `results/repeats/phase5_gene_repeat_context_summary.tsv`",
            "- `results/repeats/phase5_candidate_locus_repeat_context.tsv`",
            "- `results/evidence/phase5_repeat_context_evidence.tsv`",
            "- `results/evidence/integrated_evidence.tsv`",
        ]
    )
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_phase5_integration_outputs(
    base_evidence: Path,
    phase4e_evidence: Path,
    phase5_gene_summary: Path,
    phase5_evidence_output: Path,
    integrated_output: Path,
    report_output: Path | None = None,
    phase5c_gene_qc: Path | None = None,
) -> None:
    phase5_rows, integrated_rows = integrate_phase5_repeat_context(base_evidence, phase4e_evidence, phase5_gene_summary, phase5c_gene_qc)
    write_tsv(phase5_evidence_output, phase5_rows, EVIDENCE_COLUMNS)
    write_tsv(integrated_output, integrated_rows, EVIDENCE_COLUMNS)
    if report_output is not None:
        phase5c_rows = read_tsv(phase5c_gene_qc, PHASE5C_GENE_QC_COLUMNS) if phase5c_gene_qc and phase5c_gene_qc.exists() else None
        write_phase5_report(phase5_rows, read_tsv(phase5_gene_summary, PHASE5_GENE_REPEAT_CONTEXT_COLUMNS), phase5c_rows, report_output)
    LOGGER.info("Integrated Phase 5 repeat context for %d genes", len(phase5_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Integrate Phase 5 repeat context as evidence artifact risk.")
    parser.add_argument("--base-evidence", type=Path, required=True)
    parser.add_argument("--phase4e-evidence", type=Path, required=True)
    parser.add_argument("--phase5-gene-summary", type=Path, required=True)
    parser.add_argument("--phase5-evidence-output", type=Path, required=True)
    parser.add_argument("--integrated-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path)
    parser.add_argument("--phase5c-gene-qc", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    write_phase5_integration_outputs(
        args.base_evidence,
        args.phase4e_evidence,
        args.phase5_gene_summary,
        args.phase5_evidence_output,
        args.integrated_output,
        args.report_output,
        args.phase5c_gene_qc,
    )


if __name__ == "__main__":
    main()
