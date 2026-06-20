"""Phase 7a RNA-seq metadata and expression-readiness layer."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .schemas import (
    PHASE7_CANDIDATE_EXPRESSION_PLAN_COLUMNS,
    PHASE7_RNASEQ_MANIFEST_COLUMNS,
    PHASE7_RNASEQ_READINESS_COLUMNS,
)
from .utils import NOT_ASSESSED, clean_text, ensure_parent, read_yaml, write_tsv

LOGGER = logging.getLogger(__name__)


def get_field(row: dict, *names: str, default: str = NOT_ASSESSED) -> str:
    for name in names:
        if name in row and clean_text(row[name]):
            return clean_text(row[name])
    return default


def infer_tissue(row: dict) -> str:
    explicit = get_field(row, "tissue", "Body_Site", "body_site", default="")
    if explicit:
        return explicit
    searchable = " ".join(
        get_field(row, key, default="")
        for key in ["sample_name", "SampleName", "LibraryName", "library_name", "source"]
    ).lower()
    if "retina" in searchable:
        return "retina"
    return NOT_ASSESSED


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def parse_int(value: str) -> int:
    if value in {"", NOT_ASSESSED}:
        return 0
    return int(float(value))


def run_is_expression_ready(row: dict) -> tuple[str, bool, str, str, str]:
    strategy = get_field(row, "library_strategy", "LibraryStrategy").lower()
    source = get_field(row, "library_source", "LibrarySource").lower()
    selection = get_field(row, "library_selection", "LibrarySelection").lower()
    layout = get_field(row, "library_layout", "LibraryLayout").lower()
    tissue = infer_tissue(row).lower()
    required = [
        get_field(row, "run", "Run"),
        get_field(row, "experiment", "Experiment"),
        get_field(row, "biosample", "BioSample"),
        get_field(row, "sample_name", "SampleName"),
    ]
    if any(value in {"", NOT_ASSESSED} for value in required):
        return (
            "METADATA_INCOMPLETE_REQUIRES_REVIEW",
            False,
            "missing_run_experiment_biosample_or_sample_name",
            "high_metadata_uncertainty",
            "Run-level metadata are incomplete, so this record is not ready for expression planning.",
        )
    if strategy != "rna_seq" and strategy != "rna-seq":
        return (
            "EXCLUDED_NON_TRANSCRIPTOMIC_RESOURCE",
            False,
            "library_strategy_is_not_rna_seq",
            "not_applicable_for_expression",
            "This run is retained for provenance but excluded from candidate expression planning.",
        )
    if source != "transcriptomic":
        return (
            "EXCLUDED_NON_TRANSCRIPTOMIC_RESOURCE",
            False,
            "library_source_is_not_transcriptomic",
            "not_applicable_for_expression",
            "This run is retained for provenance but excluded from candidate expression planning.",
        )
    if tissue in {"", NOT_ASSESSED.lower()}:
        return (
            "RNA_SEQ_METADATA_REQUIRES_TISSUE_REVIEW",
            False,
            "tissue_not_resolved_from_metadata",
            "moderate_metadata_uncertainty",
            "This RNA-seq run needs tissue metadata review before candidate expression planning.",
        )
    if layout != "paired":
        return (
            "RNA_SEQ_METADATA_READY_WITH_LAYOUT_CAVEAT",
            True,
            "library_layout_is_not_paired",
            "moderate_design_risk",
            "This RNA-seq run is eligible for future candidate-gene expression audit with layout-specific handling.",
        )
    selection_note = "polyA" if selection == "polya" else selection or NOT_ASSESSED
    return (
        "RNA_SEQ_METADATA_READY_FOR_CANDIDATE_EXPRESSION_AUDIT",
        True,
        "NOT_APPLICABLE",
        "moderate_retina_only_design_risk",
        f"This {selection_note} RNA-seq run is eligible for future candidate-gene expression audit in retina; Phase 7a does not quantify expression.",
    )


def build_manifest_rows(config: dict) -> list[dict[str, str]]:
    phase7 = config.get("phase7_rnaseq", {})
    resource_id = get_field(phase7, "resource_id", default="SMIC_RETINA_PRJNA1246101_2026")
    bioproject = get_field(phase7, "bioproject", default="PRJNA1246101")
    source_url = get_field(phase7, "source_url", "runinfo_url", default="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1246101")
    retrieval_date = get_field(phase7, "retrieval_date", default=NOT_ASSESSED)
    local_root = Path(get_field(phase7, "raw_read_local_dir", default="data/raw/rnaseq/SMIC_RETINA_PRJNA1246101_2026"))
    rows: list[dict[str, str]] = []
    for run in phase7.get("registered_runs", []):
        tissue = infer_tissue(run)
        run_accession = get_field(run, "run", "Run")
        metadata_status = "RUNINFO_METADATA_REGISTERED"
        local_path = local_root / f"{run_accession}.sra"
        rows.append(
            {
                "resource_id": resource_id,
                "bioproject": get_field(run, "bioproject", "BioProject", default=bioproject),
                "sra_study": get_field(run, "sra_study", "SRAStudy"),
                "experiment": get_field(run, "experiment", "Experiment"),
                "run": run_accession,
                "biosample": get_field(run, "biosample", "BioSample"),
                "sample": get_field(run, "sample", "Sample"),
                "sample_name": get_field(run, "sample_name", "SampleName"),
                "scientific_name": get_field(run, "scientific_name", "ScientificName"),
                "taxon_id": get_field(run, "taxon_id", "TaxID"),
                "sex": get_field(run, "sex", "Sex"),
                "tissue": tissue,
                "library_strategy": get_field(run, "library_strategy", "LibraryStrategy"),
                "library_source": get_field(run, "library_source", "LibrarySource"),
                "library_selection": get_field(run, "library_selection", "LibrarySelection"),
                "library_layout": get_field(run, "library_layout", "LibraryLayout"),
                "platform": get_field(run, "platform", "Platform"),
                "model": get_field(run, "model", "Model"),
                "spots": get_field(run, "spots"),
                "bases": get_field(run, "bases"),
                "avg_length": get_field(run, "avgLength", "avg_length"),
                "size_mb": get_field(run, "size_MB", "size_mb"),
                "release_date": get_field(run, "ReleaseDate", "release_date"),
                "load_date": get_field(run, "LoadDate", "load_date"),
                "download_path": get_field(run, "download_path"),
                "source_url": source_url,
                "retrieval_date": retrieval_date,
                "local_path": str(local_path),
                "download_status": "RAW_READS_NOT_DOWNLOADED_PHASE7A_METADATA_ONLY",
                "metadata_status": metadata_status,
                "usage_notes": "RunInfo metadata only. Do not use as expression evidence until reads or processed counts are analyzed with a documented workflow.",
            }
        )
    if not rows:
        raise ValueError("config.phase7_rnaseq.registered_runs is empty")
    return rows


def build_readiness_rows(manifest_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in manifest_rows:
        status, included, exclusion, risk, interpretation = run_is_expression_ready(row)
        rows.append(
            {
                "resource_id": row["resource_id"],
                "run": row["run"],
                "experiment": row["experiment"],
                "biosample": row["biosample"],
                "sample_name": row["sample_name"],
                "tissue": row["tissue"],
                "library_strategy": row["library_strategy"],
                "library_source": row["library_source"],
                "library_selection": row["library_selection"],
                "library_layout": row["library_layout"],
                "platform_model": f"{row['platform']} {row['model']}".strip(),
                "spots": row["spots"],
                "bases": row["bases"],
                "size_mb": row["size_mb"],
                "metadata_status": row["metadata_status"],
                "expression_readiness_status": status,
                "included_for_candidate_expression_plan": bool_text(included),
                "exclusion_reason": exclusion,
                "artifact_or_design_risk": risk,
                "conservative_interpretation": interpretation,
                "required_validation": "Before expression support is reported, download raw reads or locate processed counts, record software versions and parameters, run read QC, choose a justified transcriptome/genome reference, quantify candidate genes, and inspect mapping/assignment ambiguity.",
                "claim_language_guardrail": "Use detected or expressed in retina only after quantification. Do not infer pathway activity, organism-wide state, causation, telomere activity, gene state, clinical relevance, or longevity mechanism from Phase 7a metadata.",
                "supporting_files": "data/metadata/rnaseq_manifest.tsv",
            }
        )
    return rows


def load_candidate_panel_rows(candidate_panel_path: Path) -> list[dict[str, str]]:
    panel = read_yaml(candidate_panel_path)
    panels = panel.get("panels", {})
    rows: list[dict[str, str]] = []
    for panel_id, panel_data in panels.items():
        mechanism = get_field(panel_data, "mechanism")
        for candidate in panel_data.get("candidates", []):
            synonyms = candidate.get("synonyms", [])
            rows.append(
                {
                    "mechanism": mechanism,
                    "gene_symbol": get_field(candidate, "gene_symbol"),
                    "synonyms": ";".join(synonyms) if synonyms else NOT_ASSESSED,
                    "panel_id": panel_id,
                }
            )
    if not rows:
        raise ValueError(f"No candidates found in {candidate_panel_path}")
    return rows


def recommended_method_for_plan(usable_count: int, usable_tissues: str) -> str:
    if usable_count == 0:
        return "NOT_ASSESSED_UNTIL_USABLE_RNASEQ_METADATA_OR_PROCESSED_COUNTS_ARE_AVAILABLE"
    return (
        "Future quantification should use a documented RNA-seq workflow: read QC, contamination/adaptor review, "
        "validated reference transcriptome or genome annotation, candidate-level count/TPM extraction, and manual "
        f"review of multi-mapping/paralog ambiguity in {usable_tissues}."
    )


def build_candidate_expression_plan(
    candidate_panel_path: Path,
    readiness_rows: list[dict[str, str]],
    min_replicates: int,
) -> list[dict[str, str]]:
    candidates = load_candidate_panel_rows(candidate_panel_path)
    usable = [row for row in readiness_rows if row["included_for_candidate_expression_plan"] == "true"]
    usable_experiments = ";".join(sorted({row["experiment"] for row in usable})) if usable else NOT_ASSESSED
    usable_biosamples = ";".join(sorted({row["biosample"] for row in usable})) if usable else NOT_ASSESSED
    usable_tissues = ";".join(sorted({row["tissue"] for row in usable})) if usable else NOT_ASSESSED
    resource_ids = ";".join(sorted({row["resource_id"] for row in readiness_rows})) if readiness_rows else NOT_ASSESSED
    replicate_note = (
        f"{len(usable)} usable RNA-seq runs meet the Phase 7a metadata gate."
        if len(usable) >= min_replicates
        else f"{len(usable)} usable RNA-seq runs are below the configured metadata gate of {min_replicates}."
    )
    rows: list[dict[str, str]] = []
    for candidate in candidates:
        rows.append(
            {
                "mechanism": candidate["mechanism"],
                "gene_symbol": candidate["gene_symbol"],
                "synonyms": candidate["synonyms"],
                "rnaseq_resource_id": resource_ids,
                "usable_rnaseq_run_count": str(len(usable)),
                "usable_experiments": usable_experiments,
                "usable_biosamples": usable_biosamples,
                "usable_tissues": usable_tissues,
                "candidate_panel_status": "CURATED_CANDIDATE_PANEL_ENTRY",
                "planned_quantification_scope": "candidate_gene_detection_and_expression_support_in_retina_only",
                "recommended_future_method": recommended_method_for_plan(len(usable), usable_tissues),
                "current_expression_status": "NOT_QUANTIFIED_PHASE7A_METADATA_ONLY",
                "design_limitations": f"{replicate_note} Retina-only data cannot support organism-wide aging interpretation or differential expression without a defensible design.",
                "conservative_interpretation": "This candidate is eligible for future retina RNA-seq expression audit; Phase 7a records metadata readiness only.",
                "required_inputs": "raw_FASTQ_or_processed_counts;sample_metadata;reference_transcriptome_or_genome_annotation;gene_id_mapping;read_qc",
                "required_validation": "Confirm sample identity, tissue labels, library design, replicate structure, mapping/quantification quality, and candidate gene ID mapping before reporting detected expression.",
                "claim_language_guardrail": "Do not use activated or pathway-state language. Do not generalize retina findings to the whole organism. Do not infer causation, functional advantage, clinical relevance, or longevity mechanism from Phase 7a.",
                "supporting_files": "data/metadata/rnaseq_manifest.tsv;results/rnaseq/phase7_rnaseq_readiness.tsv;config/candidate_panels.yaml",
            }
        )
    return rows


def write_report(
    manifest_rows: list[dict[str, str]],
    readiness_rows: list[dict[str, str]],
    plan_rows: list[dict[str, str]],
    output: Path,
) -> None:
    usable = [row for row in readiness_rows if row["included_for_candidate_expression_plan"] == "true"]
    excluded = [row for row in readiness_rows if row["included_for_candidate_expression_plan"] != "true"]
    total_bases = sum(parse_int(row["bases"]) for row in usable)
    lines = [
        "# Phase 7a RNA-seq Readiness Report",
        "",
        "Phase 7a records RNA-seq metadata and candidate-expression readiness only. It does not download raw reads, quantify expression, test differential expression, or report biological expression support.",
        "",
        "## Run-Level Metadata",
        "",
        f"- Registered SRA runs: {len(manifest_rows)}",
        f"- Runs eligible for future candidate expression audit: {len(usable)}",
        f"- Runs excluded from expression planning: {len(excluded)}",
        f"- Eligible tissues: {';'.join(sorted({row['tissue'] for row in usable})) if usable else NOT_ASSESSED}",
        f"- Eligible total bases: {total_bases}",
        "",
        "## Candidate Plan",
        "",
        f"- Candidate panel entries staged for future expression audit: {len(plan_rows)}",
        "- Current expression status for all candidates: `NOT_QUANTIFIED_PHASE7A_METADATA_ONLY`",
        "",
        "Guardrail: future expression results must be described as tissue-specific detection/expression support only, unless a later workflow adds a justified design and statistical model.",
        "",
        "Supporting tables:",
        "",
        "- `data/metadata/rnaseq_manifest.tsv`",
        "- `results/rnaseq/phase7_rnaseq_readiness.tsv`",
        "- `results/rnaseq/phase7_candidate_expression_plan.tsv`",
    ]
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_phase7_outputs(
    config_path: Path,
    manifest_output: Path,
    readiness_output: Path,
    expression_plan_output: Path,
    report_output: Path,
) -> None:
    config = read_yaml(config_path)
    phase7 = config.get("phase7_rnaseq", {})
    candidate_panel = Path(get_field(phase7, "candidate_panel", default="config/candidate_panels.yaml"))
    min_replicates = parse_int(get_field(phase7, "minimum_replicates_for_expression_audit", default="3"))

    manifest_rows = build_manifest_rows(config)
    readiness_rows = build_readiness_rows(manifest_rows)
    plan_rows = build_candidate_expression_plan(candidate_panel, readiness_rows, min_replicates)

    write_tsv(manifest_output, manifest_rows, PHASE7_RNASEQ_MANIFEST_COLUMNS)
    write_tsv(readiness_output, readiness_rows, PHASE7_RNASEQ_READINESS_COLUMNS)
    write_tsv(expression_plan_output, plan_rows, PHASE7_CANDIDATE_EXPRESSION_PLAN_COLUMNS)
    write_report(manifest_rows, readiness_rows, plan_rows, report_output)
    LOGGER.info(
        "Wrote Phase 7a RNA-seq readiness outputs for %d runs and %d candidate genes",
        len(manifest_rows),
        len(plan_rows),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 7a RNA-seq metadata and expression-readiness layer.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    parser.add_argument("--readiness-output", type=Path, required=True)
    parser.add_argument("--expression-plan-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    write_phase7_outputs(
        args.config,
        args.manifest_output,
        args.readiness_output,
        args.expression_plan_output,
        args.report_output,
    )


if __name__ == "__main__":
    main()
