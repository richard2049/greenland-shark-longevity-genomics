"""Phase 7b candidate expression quantification strategy."""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

from .schemas import (
    CANDIDATE_ISOFORM_AUDIT_COLUMNS,
    COPY_NUMBER_COLUMNS,
    PHASE4E_GENE_HARDENING_COLUMNS,
    PHASE4E_LOCUS_HARDENING_COLUMNS,
    PHASE7_CANDIDATE_EXPRESSION_PLAN_COLUMNS,
    PHASE7_RNASEQ_READINESS_COLUMNS,
    PHASE7B_CANDIDATE_QUANT_MAP_COLUMNS,
    PHASE7B_REFERENCE_STRATEGY_COLUMNS,
    PHASE7B_RUN_PLAN_COLUMNS,
)
from .utils import NOT_ASSESSED, ensure_parent, join_values, read_tsv, read_yaml, split_delimited_values, write_tsv

LOGGER = logging.getLogger(__name__)
HEADER_TAG_RE = re.compile(r"\[([^=\]]+)=([^\]]*)\]")


def locus_tag_from_gene_id(gene_id: str) -> str:
    if gene_id.startswith("gene-"):
        return gene_id.removeprefix("gene-")
    return gene_id


def compact_protein_id(protein_id: str) -> str:
    if not protein_id or protein_id == NOT_ASSESSED:
        return NOT_ASSESSED
    return protein_id.replace("gnl|WGS:ZZZZ|", "").replace("gnl|WGS_ZZZZ|", "")


def parse_cds_headers(cds_fasta: Path) -> list[dict[str, str]]:
    if not cds_fasta.exists():
        return []
    rows: list[dict[str, str]] = []
    with cds_fasta.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            sequence_id = header.split()[0]
            tags = {key: value for key, value in HEADER_TAG_RE.findall(header)}
            rows.append(
                {
                    "sequence_id": sequence_id,
                    "gene_symbol": tags.get("gene", NOT_ASSESSED),
                    "locus_tag": tags.get("locus_tag", NOT_ASSESSED),
                    "protein_id": tags.get("protein_id", NOT_ASSESSED),
                    "product": tags.get("protein", NOT_ASSESSED),
                }
            )
    return rows


def index_by_gene(rows: list[dict[str, str]], gene_column: str = "gene_symbol") -> dict[str, dict[str, str]]:
    return {row[gene_column]: row for row in rows}


def group_isoforms(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    grouped: dict[str, dict[str, str]] = {}
    for row in rows:
        grouped[row["gene_symbol"]] = row
    return grouped


def collect_phase4e_candidate_gene_ids(locus_rows: list[dict[str, str]]) -> dict[str, set[str]]:
    by_gene: dict[str, set[str]] = {}
    for row in locus_rows:
        gene = row["gene_symbol"]
        by_gene.setdefault(gene, set()).update(split_delimited_values(row.get("overlapping_gene_ids", "")))
    return by_gene


def cds_by_locus_tag(cds_rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    by_locus: dict[str, list[dict[str, str]]] = {}
    for row in cds_rows:
        locus = row.get("locus_tag", NOT_ASSESSED)
        if locus and locus != NOT_ASSESSED:
            by_locus.setdefault(locus, []).append(row)
    return by_locus


def local_input_status(paths: list[Path]) -> str:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        return "MISSING:" + ";".join(missing)
    return "ALL_REQUIRED_LOCAL_INPUTS_PRESENT"


def build_reference_strategy_rows(config: dict) -> list[dict[str, str]]:
    phase7b = config["phase7b_expression_quantification_strategy"]
    genome = Path(phase7b["genome_fasta"])
    gff = Path(phase7b["annotation_gff"])
    cds = Path(phase7b["cds_fasta"])
    transcriptome = Path(phase7b["transcriptome_fasta_planned"])
    tx2gene = Path(phase7b["tx2gene_planned"])
    salmon_index = Path(phase7b["salmon_index_planned"])
    candidate_matrix = Path(phase7b["candidate_quantification_matrix_planned"])
    return [
        {
            "strategy_id": "PHASE7B_PRIMARY_SALMON_SELECTIVE_ALIGNMENT",
            "strategy_role": "preferred_future_quantification",
            "recommended_tool_or_method": "Salmon with validateMappings/selective-alignment-style indexing and transcript-to-gene summarization",
            "algorithm_class": "alignment-aware transcript-level quantification with probabilistic multi-mapping assignment",
            "required_inputs": f"{transcriptome};{tx2gene};{genome};raw paired-end FASTQ files",
            "local_inputs_available": local_input_status([genome, gff]),
            "planned_output": str(candidate_matrix),
            "readiness_status": "PLANNED_REQUIRES_RAW_READS_AND_TRANSCRIPT_FASTA_GENERATION",
            "rationale": "Efficient for paired-end RNA-seq and appropriate for candidate-level TPM/count extraction while retaining transcript-to-gene ambiguity information.",
            "limitations": "Requires generated transcript FASTA, tx2gene map, raw-read QC, and careful handling of duplicated or paralogous candidates.",
            "claim_language_guardrail": "Future results may support detected/expression in retina only; do not infer activation, pathway state, causation, or longevity mechanism.",
            "next_action": "Generate transcript FASTA and tx2gene from the Tokyo genome/GFF, build decoy-aware Salmon index, then quantify each metadata-ready RNA-seq run.",
        },
        {
            "strategy_id": "PHASE7B_LOCAL_CDS_ONLY_PILOT",
            "strategy_role": "coding_sequence_only_fallback",
            "recommended_tool_or_method": "Candidate tx2gene mapping from complete.cds.fna headers; optional pilot Salmon index on CDS sequences",
            "algorithm_class": "coding-sequence reference mapping",
            "required_inputs": str(cds),
            "local_inputs_available": local_input_status([cds]),
            "planned_output": "results/rnaseq/phase7b_candidate_quantification_map.tsv",
            "readiness_status": "LOCAL_CDS_REFERENCE_AVAILABLE_CODING_ONLY",
            "rationale": "The CDS FASTA is local, small, and linked to gene/locus/protein identifiers, so it can validate candidate-to-reference mapping before raw-read processing.",
            "limitations": "CDS-only quantification would ignore UTR reads and should not replace a transcriptome or genome-aware final expression audit.",
            "claim_language_guardrail": "Do not treat CDS-only readiness as expression evidence.",
            "next_action": "Use CDS headers to build candidate reference mapping and identify candidates needing manual transcript/locus review.",
        },
        {
            "strategy_id": "PHASE7B_GENOME_ALIGNED_VALIDATION",
            "strategy_role": "validation_for_ambiguous_candidates",
            "recommended_tool_or_method": "STAR or HISAT2 plus featureCounts or equivalent locus-aware counting",
            "algorithm_class": "splice-aware genome alignment and feature-level counting",
            "required_inputs": f"{genome};{gff};raw paired-end FASTQ files",
            "local_inputs_available": local_input_status([genome, gff]),
            "planned_output": "results/rnaseq/phase7b_high_ambiguity_locus_validation.tsv",
            "readiness_status": "PLANNED_FOR_PARALOG_OR_REPEAT_SENSITIVE_CANDIDATES",
            "rationale": "Genome-aligned review is better for checking multi-mapping, locus specificity, and gene-family ambiguity around FTH1B, H1F0, RAD51, and TP53.",
            "limitations": "Heavier than Salmon and still requires manual review of mapping quality, multi-mapping policy, and annotation overlap.",
            "claim_language_guardrail": "Genome-aligned counts still do not prove pathway activity, functional advantage, or organism-wide longevity relevance.",
            "next_action": "Apply only to candidates with high ambiguity or after Salmon flags multi-mapping/paralog uncertainty.",
        },
    ]


def risk_for_candidate(
    gene: str,
    copy_row: dict[str, str] | None,
    isoform_row: dict[str, str] | None,
    phase4e_row: dict[str, str] | None,
    high_ambiguity_genes: set[str],
) -> str:
    if gene in high_ambiguity_genes:
        return "high_candidate_family_or_locus_ambiguity"
    if phase4e_row and phase4e_row.get("artifact_risk") == "high":
        return "high_phase4e_artifact_risk"
    if copy_row and copy_row.get("mapping_status", "").startswith("ANNOTATION_UNCERTAINTY"):
        return "high_annotation_uncertainty"
    if copy_row and int(copy_row.get("copy_count", "0") or 0) > 1:
        return "moderate_copy_number_or_paralog_risk"
    if isoform_row and int(isoform_row.get("isoform_count", "0") or 0) > 1:
        return "moderate_isoform_summarization_risk"
    return "low_first_pass_quantification_risk"


def status_for_candidate(
    gene_ids: list[str],
    cds_ids: list[str],
    risk: str,
    copy_row: dict[str, str] | None,
    phase4e_row: dict[str, str] | None,
) -> str:
    if not gene_ids:
        return "NOT_READY_NO_CANDIDATE_GENE_MODEL_MAP"
    if not cds_ids:
        return "NOT_READY_NO_CDS_REFERENCE_SEQUENCE_FOR_MAPPED_GENE"
    if phase4e_row and phase4e_row.get("evidence_tier_recommendation") == "Artifact/uncertain":
        return "NOT_READY_TARGETED_LOCUS_VALIDATION_REQUIRED_BEFORE_EXPRESSION_INTERPRETATION"
    if risk.startswith("high"):
        return "TARGETED_LOCUS_VALIDATION_REQUIRED_AFTER_FIRST_PASS_QUANTIFICATION"
    if risk.startswith("moderate"):
        return "READY_WITH_ISOFORM_OR_PARALOG_SUMMARIZATION_CAVEAT"
    if copy_row and copy_row.get("mapping_status", "").startswith("ANNOTATION_SYMBOL_MATCH"):
        return "READY_FOR_FIRST_PASS_CANDIDATE_QUANTIFICATION"
    return "READY_WITH_MAPPING_CAVEAT"


def scope_for_status(status: str) -> str:
    if status.startswith("READY_FOR_FIRST_PASS"):
        return "candidate_gene_detection_and_expression_support_in_retina_using_transcript_to_gene_summarization"
    if status.startswith("READY_WITH"):
        return "candidate_gene_detection_in_retina_with_isoform_or_paralog_caveats"
    if status.startswith("TARGETED"):
        return "candidate_locus_screening_only_until_manual_locus_validation_resolves_ambiguity"
    return "do_not_quantify_as_candidate_expression_support_until_reference_mapping_is_resolved"


def build_candidate_quant_map(config: dict) -> list[dict[str, str]]:
    phase7b = config["phase7b_expression_quantification_strategy"]
    plan_rows = read_tsv(Path(phase7b["phase7a_candidate_expression_plan"]), PHASE7_CANDIDATE_EXPRESSION_PLAN_COLUMNS)
    copy_rows = read_tsv(Path(phase7b["candidate_copy_number"]), COPY_NUMBER_COLUMNS)
    isoform_rows = read_tsv(Path(phase7b["candidate_isoform_audit"]), CANDIDATE_ISOFORM_AUDIT_COLUMNS)
    phase4e_gene_rows = read_tsv(Path(phase7b["phase4e_gene_summary"]), PHASE4E_GENE_HARDENING_COLUMNS)
    phase4e_locus_rows = read_tsv(Path(phase7b["phase4e_locus_review"]), PHASE4E_LOCUS_HARDENING_COLUMNS)
    cds_rows = parse_cds_headers(Path(phase7b["cds_fasta"]))

    copy_by_gene = index_by_gene(copy_rows)
    isoform_by_gene = group_isoforms(isoform_rows)
    phase4e_by_gene = index_by_gene(phase4e_gene_rows)
    phase4e_gene_ids = collect_phase4e_candidate_gene_ids(phase4e_locus_rows)
    cds_lookup = cds_by_locus_tag(cds_rows)
    high_ambiguity = set(phase7b.get("high_ambiguity_genes", []))

    rows: list[dict[str, str]] = []
    for plan in plan_rows:
        gene = plan["gene_symbol"]
        copy_row = copy_by_gene.get(gene)
        isoform_row = isoform_by_gene.get(gene)
        phase4e_row = phase4e_by_gene.get(gene)
        exact_gene_ids = split_delimited_values(copy_row.get("gene_ids", "") if copy_row else "")
        rescue_gene_ids = sorted(phase4e_gene_ids.get(gene, set()))
        quant_gene_ids = exact_gene_ids if exact_gene_ids else rescue_gene_ids
        locus_tags = [locus_tag_from_gene_id(gene_id) for gene_id in quant_gene_ids]
        matched_cds = []
        for locus_tag in locus_tags:
            matched_cds.extend(cds_lookup.get(locus_tag, []))
        cds_ids = sorted({row["sequence_id"] for row in matched_cds})
        compact_proteins = sorted({compact_protein_id(row["protein_id"]) for row in matched_cds if row["protein_id"]})
        if isoform_row and isoform_row.get("representative_protein_id"):
            compact_proteins.append(compact_protein_id(isoform_row["representative_protein_id"]))
        compact_proteins = sorted({value for value in compact_proteins if value and value != NOT_ASSESSED})

        risk = risk_for_candidate(gene, copy_row, isoform_row, phase4e_row, high_ambiguity)
        status = status_for_candidate(quant_gene_ids, cds_ids, risk, copy_row, phase4e_row)
        rows.append(
            {
                "mechanism": plan["mechanism"],
                "gene_symbol": gene,
                "phase7a_current_expression_status": plan["current_expression_status"],
                "phase7a_usable_run_count": plan["usable_rnaseq_run_count"],
                "phase7a_usable_tissues": plan["usable_tissues"],
                "exact_annotation_gene_ids": join_values(exact_gene_ids),
                "phase4e_candidate_gene_ids": join_values(rescue_gene_ids),
                "quantification_gene_ids": join_values(quant_gene_ids),
                "cds_reference_sequence_ids": join_values(cds_ids),
                "transcript_or_cds_count": str(len(cds_ids)),
                "representative_protein_ids": join_values(compact_proteins),
                "isoform_count": isoform_row.get("isoform_count", NOT_ASSESSED) if isoform_row else NOT_ASSESSED,
                "orthology_mapping_status": copy_row.get("mapping_status", NOT_ASSESSED) if copy_row else NOT_ASSESSED,
                "phase4e_status": phase4e_row.get("phase4e_hardened_status", NOT_ASSESSED) if phase4e_row else NOT_ASSESSED,
                "artifact_or_ambiguity_risk": risk,
                "quantification_strategy_status": status,
                "recommended_quantification_scope": scope_for_status(status),
                "required_validation": "Download/validate reads, run read QC, generate transcript reference and tx2gene map, quantify with recorded parameters, inspect multi-mapping/paralog ambiguity, and report retina-specific detected/expression status only.",
                "claim_language_guardrail": "Do not infer activation, differential expression, pathway activity, causation, functional advantage, clinical relevance, or longevity mechanism from Phase 7b strategy outputs.",
                "supporting_files": "results/rnaseq/phase7_candidate_expression_plan.tsv;results/orthology/candidate_copy_number.tsv;results/validation/candidate_isoform_audit.tsv;results/rescue/phase4e_gene_hardened_summary.tsv;complete.cds.fna",
            }
        )
    return rows


def build_run_plan_rows(config: dict) -> list[dict[str, str]]:
    phase7b = config["phase7b_expression_quantification_strategy"]
    runs = [run["run"] for run in config["phase7_rnaseq"]["registered_runs"] if run.get("library_strategy") == "RNA-Seq"]
    run_list = ",".join(runs)
    return [
        {
            "step_id": "01",
            "step_name": "download_or_link_raw_reads",
            "recommended_environment": "WSL/Linux or Docker; avoid broad native-Windows Snakemake",
            "command_template": f"prefetch {run_list} && fasterq-dump --split-files --threads <threads> <SRR>",
            "expected_output": "data/raw/rnaseq/SMIC_RETINA_PRJNA1246101_2026/<SRR>_1.fastq.gz and <SRR>_2.fastq.gz",
            "run_status": "PLANNED_NOT_RUN",
            "notes": "Raw reads are not downloaded in Phase 7b. Record checksums, SRA Toolkit version, and read counts if this step is run.",
        },
        {
            "step_id": "02",
            "step_name": "read_qc",
            "recommended_environment": "WSL/Linux or Docker",
            "command_template": "fastp or FastQC/MultiQC on each paired-end run with logs retained",
            "expected_output": "results/rnaseq/qc/<SRR> read-quality reports",
            "run_status": "PLANNED_NOT_RUN",
            "notes": "QC must precede expression interpretation; low-quality or contaminated reads should block expression support.",
        },
        {
            "step_id": "03",
            "step_name": "generate_transcript_reference",
            "recommended_environment": "WSL/Linux or Docker",
            "command_template": f"gffread {phase7b['annotation_gff']} -g {phase7b['genome_fasta']} -w {phase7b['transcriptome_fasta_planned']}",
            "expected_output": f"{phase7b['transcriptome_fasta_planned']};{phase7b['tx2gene_planned']}",
            "run_status": "PLANNED_NOT_RUN",
            "notes": "Prefer full transcript sequences from genome+GFF. complete.cds.fna remains a coding-only fallback.",
        },
        {
            "step_id": "04",
            "step_name": "build_salmon_index",
            "recommended_environment": "WSL/Linux or Docker",
            "command_template": "salmon index --validateMappings-style selective-alignment/decoy-aware setup using generated transcript FASTA and genome decoys",
            "expected_output": str(phase7b["salmon_index_planned"]),
            "run_status": "PLANNED_NOT_RUN",
            "notes": "Record Salmon version, index parameters, transcript FASTA checksum, and decoy policy.",
        },
        {
            "step_id": "05",
            "step_name": "quantify_each_run",
            "recommended_environment": "WSL/Linux or Docker",
            "command_template": "salmon quant -i <index> -l A -1 <SRR>_1.fastq.gz -2 <SRR>_2.fastq.gz --validateMappings --seqBias --gcBias",
            "expected_output": "results/rnaseq/salmon/<SRR>/quant.sf",
            "run_status": "PLANNED_NOT_RUN",
            "notes": "This supports candidate detection/expression summaries, not differential expression by itself.",
        },
        {
            "step_id": "06",
            "step_name": "candidate_summarization_and_ambiguity_review",
            "recommended_environment": "green_shark Python environment",
            "command_template": "future repository module to join quant.sf, tx2gene, phase7b_candidate_quantification_map.tsv, and high-ambiguity flags",
            "expected_output": str(phase7b["candidate_quantification_matrix_planned"]),
            "run_status": "PLANNED_NOT_RUN",
            "notes": "High-risk genes require manual multi-mapping/locus review before expression support is used in Phase 8.",
        },
    ]


def write_report(
    reference_rows: list[dict[str, str]],
    candidate_rows: list[dict[str, str]],
    run_plan_rows: list[dict[str, str]],
    output: Path,
) -> None:
    ready = [row for row in candidate_rows if row["quantification_strategy_status"].startswith("READY")]
    high_risk = [row for row in candidate_rows if row["artifact_or_ambiguity_risk"].startswith("high")]
    not_ready = [row for row in candidate_rows if row["quantification_strategy_status"].startswith("NOT_READY")]
    lines = [
        "# Phase 7b Candidate Expression Quantification Strategy",
        "",
        "Phase 7b defines the quantification strategy and candidate-to-reference mapping. It does not download reads, run Salmon, quantify expression, or test differential expression.",
        "",
        "## Algorithm Choice",
        "",
        "The recommended primary route is Salmon-style transcript quantification using an annotation-derived transcript FASTA with genome decoys. This is efficient for three paired-end retinal RNA-seq runs and preserves transcript/gene assignment ambiguity. Genome-aligned counting is reserved as a validation layer for duplicated, paralogous, repeat-overlapping, or unresolved loci.",
        "",
        "## Candidate Readiness",
        "",
        f"- Candidate-panel genes evaluated: {len(candidate_rows)}",
        f"- Ready or ready with caveat for first-pass quantification: {len(ready)}",
        f"- High-ambiguity candidates requiring targeted review: {len(high_risk)}",
        f"- Not ready for expression-support interpretation: {len(not_ready)}",
        "",
        "## Planned Steps",
        "",
    ]
    for row in run_plan_rows:
        lines.append(f"- {row['step_id']}. {row['step_name']}: {row['run_status']}")
    lines.extend(
        [
            "",
            "Guardrail: future results may be reported only as retina-specific detected/expression support after read QC, quantification, gene-ID mapping, and ambiguity review.",
            "",
            "Supporting tables:",
            "",
            "- `results/rnaseq/phase7b_reference_quantification_strategy.tsv`",
            "- `results/rnaseq/phase7b_candidate_quantification_map.tsv`",
            "- `results/rnaseq/phase7b_quantification_run_plan.tsv`",
        ]
    )
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_phase7b_outputs(
    config_path: Path,
    reference_strategy_output: Path,
    candidate_quant_map_output: Path,
    run_plan_output: Path,
    report_output: Path,
) -> None:
    config = read_yaml(config_path)
    if "phase7b_expression_quantification_strategy" not in config:
        raise ValueError("config.yaml is missing phase7b_expression_quantification_strategy")
    reference_rows = build_reference_strategy_rows(config)
    candidate_rows = build_candidate_quant_map(config)
    run_plan_rows = build_run_plan_rows(config)
    write_tsv(reference_strategy_output, reference_rows, PHASE7B_REFERENCE_STRATEGY_COLUMNS)
    write_tsv(candidate_quant_map_output, candidate_rows, PHASE7B_CANDIDATE_QUANT_MAP_COLUMNS)
    write_tsv(run_plan_output, run_plan_rows, PHASE7B_RUN_PLAN_COLUMNS)
    write_report(reference_rows, candidate_rows, run_plan_rows, report_output)
    LOGGER.info("Wrote Phase 7b strategy outputs for %d candidates", len(candidate_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 7b candidate expression quantification strategy.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--reference-strategy-output", type=Path, required=True)
    parser.add_argument("--candidate-quant-map-output", type=Path, required=True)
    parser.add_argument("--run-plan-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    write_phase7b_outputs(
        args.config,
        args.reference_strategy_output,
        args.candidate_quant_map_output,
        args.run_plan_output,
        args.report_output,
    )


if __name__ == "__main__":
    main()
