"""Phase 7e hardening of retina RNA-seq expression interpretation."""

from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path
from statistics import median

from .schemas import (
    DUPLICATION_AUDIT_COLUMNS,
    PHASE4E_GENE_HARDENING_COLUMNS,
    PHASE5C_GENE_QC_COLUMNS,
    PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS,
    PHASE7D_CANDIDATE_EXPRESSION_MATRIX_COLUMNS,
    PHASE7D_FASTQ_QC_COLUMNS,
    PHASE7D_SALMON_QUANT_SUMMARY_COLUMNS,
    PHASE7E_CANDIDATE_HARDENED_COLUMNS,
    PHASE7E_PARAMETER_REVIEW_COLUMNS,
    PHASE7E_RUN_QC_COLUMNS,
)
from .utils import NOT_ASSESSED, ensure_parent, join_values, read_tsv, read_yaml, write_tsv

LOGGER = logging.getLogger(__name__)
DECOY_RE = re.compile(r"Number of decoys\s*:\s*(\d+)")


def parse_float(value: str | int | float | None) -> float | None:
    if value in {None, "", NOT_ASSESSED}:
        return None
    return float(str(value).replace(",", "."))


def parse_int(value: str | int | float | None) -> int:
    parsed = parse_float(value)
    return 0 if parsed is None else int(parsed)


def format_float(value: float | None, digits: int = 6) -> str:
    if value is None:
        return NOT_ASSESSED
    return f"{value:.{digits}f}"


def index_by_gene(rows: list[dict[str, str]], column: str = "gene_symbol") -> dict[str, dict[str, str]]:
    return {row[column]: row for row in rows}


def group_by(rows: list[dict[str, str]], column: str) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row[column], []).append(row)
    return grouped


def load_salmon_meta(output_dir: Path) -> dict[str, str]:
    meta_path = output_dir / "aux_info" / "meta_info.json"
    if not meta_path.exists():
        return {
            "num_processed": NOT_ASSESSED,
            "num_mapped": NOT_ASSESSED,
            "percent_mapped": NOT_ASSESSED,
            "library_types": NOT_ASSESSED,
            "meta_info": str(meta_path),
        }
    with meta_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    library_types = data.get("library_types", [])
    return {
        "num_processed": str(data.get("num_processed", NOT_ASSESSED)),
        "num_mapped": str(data.get("num_mapped", NOT_ASSESSED)),
        "percent_mapped": str(data.get("percent_mapped", NOT_ASSESSED)),
        "library_types": ";".join(library_types) if library_types else NOT_ASSESSED,
        "meta_info": str(meta_path),
    }


def parse_decoy_count(stderr_log: Path) -> str:
    if not stderr_log.exists():
        return NOT_ASSESSED
    text = stderr_log.read_text(encoding="utf-8", errors="replace")
    match = DECOY_RE.search(text)
    return match.group(1) if match else NOT_ASSESSED


def mapping_rate_status(percent_mapped: float | None, low_warning: float, very_low: float) -> str:
    if percent_mapped is None:
        return "MAPPING_RATE_NOT_ASSESSED"
    if percent_mapped < very_low:
        return "VERY_LOW_MAPPING_RATE_REVIEW_REQUIRED"
    if percent_mapped < low_warning:
        return "LOW_MAPPING_RATE_CAUTIOUS_USE"
    return "PASS_FIRST_PASS_MAPPING_RATE"


def run_interpretation_use_status(quant_status: str, fastq_status: str, mapping_status: str) -> str:
    if quant_status != "SALMON_QUANT_COMPLETE":
        return "DO_NOT_USE_NOT_QUANTIFIED"
    if fastq_status != "PASS_LIGHTWEIGHT_FASTQ_QC":
        return "DO_NOT_USE_FASTQ_QC_FAILED"
    if mapping_status == "VERY_LOW_MAPPING_RATE_REVIEW_REQUIRED":
        return "LIMITED_USE_MAPPING_RATE_REVIEW_REQUIRED"
    if mapping_status == "LOW_MAPPING_RATE_CAUTIOUS_USE":
        return "USE_WITH_LOW_MAPPING_RATE_CAVEAT"
    return "USE_FOR_FIRST_PASS_RETINA_EXPRESSION_REVIEW"


def build_run_qc_rows(config: dict, quant_rows: list[dict[str, str]], fastq_qc_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    phase = config["phase7e_expression_hardening"]
    low_warning = float(phase["low_mapping_rate_warning_percent"])
    very_low = float(phase["very_low_mapping_rate_percent"])
    fastq_by_run = group_by(fastq_qc_rows, "run")
    rows: list[dict[str, str]] = []
    for quant in quant_rows:
        run = quant["run"]
        fastq_rows = fastq_by_run.get(run, [])
        qc_status_values = {row["qc_status"] for row in fastq_rows}
        fastq_status = "PASS_LIGHTWEIGHT_FASTQ_QC" if qc_status_values == {"PASS_LIGHTWEIGHT_FASTQ_QC"} and len(fastq_rows) == 2 else "FASTQ_QC_REVIEW_REQUIRED"
        read_lengths = [parse_float(row["mean_read_length"]) for row in fastq_rows]
        read_lengths = [value for value in read_lengths if value is not None]
        qualities = [parse_float(row["mean_phred_quality"]) for row in fastq_rows]
        qualities = [value for value in qualities if value is not None]
        n_percents = [parse_float(row["n_percent"]) for row in fastq_rows]
        n_percents = [value for value in n_percents if value is not None]
        records = join_values([row["records_inspected"] for row in fastq_rows])

        output_dir = Path(quant["salmon_output_dir"])
        meta = load_salmon_meta(output_dir)
        percent_mapped = parse_float(meta["percent_mapped"])
        map_status = mapping_rate_status(percent_mapped, low_warning, very_low)
        use_status = run_interpretation_use_status(quant["quantification_status"], fastq_status, map_status)
        stderr_log = Path(quant["log_stderr"])
        decoy_count = parse_decoy_count(stderr_log)
        required_validation = (
            "Review Salmon mapping rate, library type, index decoy status, and candidate-level ambiguity before Phase 8 use."
            if use_status != "USE_FOR_FIRST_PASS_RETINA_EXPRESSION_REVIEW"
            else "Retain retina-only interpretation and review candidate-specific ambiguity before Phase 8."
        )
        rows.append(
            {
                "run": run,
                "tissue": quant["tissue"],
                "quantification_status": quant["quantification_status"],
                "fastq_qc_status": fastq_status,
                "fastq_files_reviewed": str(len(fastq_rows)),
                "records_inspected_per_read": records,
                "mean_read_length_range": f"{min(read_lengths):.2f}-{max(read_lengths):.2f}" if read_lengths else NOT_ASSESSED,
                "min_mean_phred_quality": format_float(min(qualities), 3) if qualities else NOT_ASSESSED,
                "max_n_percent": format_float(max(n_percents), 3) if n_percents else NOT_ASSESSED,
                "salmon_num_processed": meta["num_processed"],
                "salmon_num_mapped": meta["num_mapped"],
                "salmon_percent_mapped": format_float(percent_mapped, 6),
                "salmon_library_types": meta["library_types"],
                "salmon_decoy_count": decoy_count,
                "mapping_rate_status": map_status,
                "interpretation_use_status": use_status,
                "required_validation": required_validation,
                "claim_language_guardrail": "Run-level RNA-seq QC is technical evidence only. Do not infer pathway activity, organism-wide aging state, causation, or longevity mechanism.",
                "supporting_files": join_values([quant["quant_sf"], quant["log_stderr"], meta["meta_info"]]),
            }
        )
    return rows


def expression_status_for_gene(
    phase7c_readiness: str,
    detected_count: int,
    quantified_count: int,
    locus_review_count: int,
    min_detected_runs: int,
    min_detection_fraction: float,
    phase4e_row: dict[str, str] | None,
) -> tuple[str, str]:
    if quantified_count == 0:
        return "NOT_QUANTIFIED", "DO_NOT_USE_AS_PHASE8_EXPRESSION_SUPPORT"
    if phase7c_readiness.startswith("NOT_QUANTIFIABLE"):
        return "EXPRESSION_NOT_INTERPRETABLE_REFERENCE_NOT_QUANTIFIABLE", "DO_NOT_USE_AS_PHASE8_EXPRESSION_SUPPORT"
    if "NOT_READY" in phase7c_readiness or "TARGETED_LOCUS" in phase7c_readiness or locus_review_count > 0 or phase4e_row is not None:
        return "LOCUS_REVIEW_REQUIRED_BEFORE_EXPRESSION_SUPPORT", "DEFER_PHASE8_EXPRESSION_SUPPORT_UNTIL_LOCUS_REVIEW"
    detection_fraction = detected_count / quantified_count if quantified_count else 0.0
    if detected_count >= min_detected_runs and detection_fraction >= min_detection_fraction:
        return "RETINA_DETECTED_FIRST_PASS_WITH_QC_CAVEATS", "CAUTIOUS_RETINA_EXPRESSION_SUPPORT_FOR_PHASE8"
    if detected_count > 0:
        return "RETINA_DETECTED_IN_LIMITED_RUNS_EXPLORATORY", "LIMITED_EXPLORATORY_EXPRESSION_SUPPORT"
    return "NOT_DETECTED_UNDER_CURRENT_RETINA_QUANTIFICATION", "NO_RETINA_EXPRESSION_SUPPORT_UNDER_CURRENT_REFERENCE"


def combine_artifact_risk(*values: str) -> str:
    parts: list[str] = []
    for value in values:
        if value in {"", NOT_ASSESSED, None}:
            continue
        for part in str(value).split(";"):
            cleaned = part.strip()
            if cleaned and cleaned != NOT_ASSESSED and cleaned not in parts:
                parts.append(cleaned)
    return ";".join(parts) if parts else NOT_ASSESSED


def conservative_interpretation(
    status: str,
    gene: str,
    detected_count: int,
    quantified_count: int,
    median_tpm_value: float | None,
    low_mapping_count: int,
    very_low_mapping_count: int,
) -> str:
    tpm_text = format_float(median_tpm_value, 6)
    qc_text = f" Run-level caveats include {low_mapping_count} low-mapping run(s), including {very_low_mapping_count} very-low-mapping run(s)."
    if status == "RETINA_DETECTED_FIRST_PASS_WITH_QC_CAVEATS":
        return (
            f"{gene} is detected in retina in {detected_count}/{quantified_count} quantified runs in this first-pass Salmon matrix "
            f"(median TPM {tpm_text}).{qc_text} Treat this as tissue-specific expression support only."
        )
    if status == "RETINA_DETECTED_IN_LIMITED_RUNS_EXPLORATORY":
        return (
            f"{gene} is detected in only {detected_count}/{quantified_count} quantified retina runs under the current thresholds "
            f"(median TPM {tpm_text}).{qc_text} Treat this as exploratory and not robust expression support."
        )
    if status == "LOCUS_REVIEW_REQUIRED_BEFORE_EXPRESSION_SUPPORT":
        return (
            f"{gene} has quantified retina signal but candidate/reference or locus ambiguity prevents using it as expression support until the locus model is reviewed."
            f"{qc_text}"
        )
    if status == "EXPRESSION_NOT_INTERPRETABLE_REFERENCE_NOT_QUANTIFIABLE":
        return f"{gene} is not interpretable in Phase 7e because the candidate could not be mapped to a usable expression reference."
    if status == "NOT_DETECTED_UNDER_CURRENT_RETINA_QUANTIFICATION":
        return (
            f"{gene} is not detected under the current retina quantification thresholds and reference model. "
            "This is condition-specific and must not be interpreted as gene absence or lack of biological function."
        )
    return f"{gene} has no completed quantification available for expression interpretation."


def build_candidate_hardened_rows(
    config: dict,
    expression_rows: list[dict[str, str]],
    phase7c_rows: list[dict[str, str]],
    duplication_rows: list[dict[str, str]],
    phase4e_rows: list[dict[str, str]],
    phase5c_rows: list[dict[str, str]],
    run_qc_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    phase = config["phase7e_expression_hardening"]
    min_detected_runs = int(phase["min_detected_runs_for_consistent_retina_support"])
    min_detection_fraction = float(phase["min_detection_fraction_for_consistent_retina_support"])
    grouped = group_by(expression_rows, "gene_symbol")
    phase7c_by_gene = index_by_gene(phase7c_rows)
    duplication_by_gene = index_by_gene(duplication_rows)
    phase4e_by_gene = index_by_gene(phase4e_rows)
    phase5c_by_gene = index_by_gene(phase5c_rows)
    run_qc_by_run = {row["run"]: row for row in run_qc_rows}

    rows: list[dict[str, str]] = []
    for gene, gene_rows in grouped.items():
        first = gene_rows[0]
        phase7c = phase7c_by_gene.get(gene, {})
        duplication = duplication_by_gene.get(gene, {})
        phase4e = phase4e_by_gene.get(gene)
        phase5c = phase5c_by_gene.get(gene, {})

        quantified = [row for row in gene_rows if row["quantification_status"] == "SALMON_QUANT_COMPLETE"]
        detected = [row for row in quantified if row["candidate_expression_status"] == "DETECTED_IN_RETINA_EXPLORATORY"]
        not_detected = [row for row in quantified if row["candidate_expression_status"] == "NOT_DETECTED_UNDER_CURRENT_RETINA_QUANTIFICATION"]
        locus_review = [row for row in quantified if row["candidate_expression_status"] == "QUANTIFIED_REQUIRES_LOCUS_REVIEW_BEFORE_INTERPRETATION"]
        tpm_values = [parse_float(row["summed_tpm"]) for row in quantified]
        tpm_values = [value for value in tpm_values if value is not None]
        read_values = [parse_float(row["summed_numreads"]) for row in quantified]
        read_values = [value for value in read_values if value is not None]
        run_qc = [run_qc_by_run[row["run"]] for row in quantified if row["run"] in run_qc_by_run]
        low_mapping_count = sum(1 for row in run_qc if row["mapping_rate_status"] in {"LOW_MAPPING_RATE_CAUTIOUS_USE", "VERY_LOW_MAPPING_RATE_REVIEW_REQUIRED"})
        very_low_mapping_count = sum(1 for row in run_qc if row["mapping_rate_status"] == "VERY_LOW_MAPPING_RATE_REVIEW_REQUIRED")
        mapping_rates = [parse_float(row["salmon_percent_mapped"]) for row in run_qc]
        mapping_rates = [value for value in mapping_rates if value is not None]
        phase7c_readiness = phase7c.get("phase7c_quantification_readiness", first["phase7c_quantification_readiness"])
        status, support_level = expression_status_for_gene(
            phase7c_readiness,
            len(detected),
            len(quantified),
            len(locus_review),
            min_detected_runs,
            min_detection_fraction,
            phase4e,
        )
        artifact_risk = combine_artifact_risk(
            phase7c.get("ambiguity_status", NOT_ASSESSED),
            duplication.get("artifact_risk", NOT_ASSESSED),
            phase4e.get("artifact_risk", NOT_ASSESSED) if phase4e else NOT_ASSESSED,
            phase5c.get("artifact_risk_modifier", NOT_ASSESSED),
            "low_mapping_rate_caveat" if low_mapping_count else NOT_ASSESSED,
            "very_low_mapping_rate_caveat" if very_low_mapping_count else NOT_ASSESSED,
        )
        required_validation = combine_artifact_risk(
            phase7c.get("required_validation", NOT_ASSESSED),
            duplication.get("required_validation", NOT_ASSESSED),
            phase4e.get("required_validation", NOT_ASSESSED) if phase4e else NOT_ASSESSED,
            phase5c.get("required_validation", NOT_ASSESSED),
            "Confirm expression-support interpretation with decoy-aware Salmon or genome-aligned counting for ambiguous loci before Phase 8 claims.",
        )
        rows.append(
            {
                "mechanism": first["mechanism"],
                "gene_symbol": gene,
                "tissue": first["tissue"],
                "quantified_run_count": str(len(quantified)),
                "detected_run_count": str(len(detected)),
                "not_detected_run_count": str(len(not_detected)),
                "locus_review_run_count": str(len(locus_review)),
                "detected_runs": join_values([row["run"] for row in detected]),
                "median_tpm": format_float(median(tpm_values) if tpm_values else None, 6),
                "min_tpm": format_float(min(tpm_values) if tpm_values else None, 6),
                "max_tpm": format_float(max(tpm_values) if tpm_values else None, 6),
                "median_numreads": format_float(median(read_values) if read_values else None, 3),
                "phase7c_quantification_readiness": phase7c_readiness,
                "reference_mapping_status": phase7c.get("reference_mapping_status", NOT_ASSESSED),
                "reference_ambiguity_status": phase7c.get("ambiguity_status", NOT_ASSESSED),
                "matched_transcript_count": phase7c.get("matched_transcript_count", first["matched_transcript_count"]),
                "orthology_mapping_status": duplication.get("mapping_status", NOT_ASSESSED),
                "duplication_artifact_risk": duplication.get("artifact_risk", NOT_ASSESSED),
                "phase4e_hardened_status": phase4e.get("phase4e_hardened_status", NOT_ASSESSED) if phase4e else NOT_ASSESSED,
                "phase4e_artifact_risk": phase4e.get("artifact_risk", NOT_ASSESSED) if phase4e else NOT_ASSESSED,
                "phase5c_artifact_risk_modifier": phase5c.get("artifact_risk_modifier", NOT_ASSESSED),
                "low_mapping_run_count": str(low_mapping_count),
                "very_low_mapping_run_count": str(very_low_mapping_count),
                "min_run_mapping_rate": format_float(min(mapping_rates) if mapping_rates else None, 6),
                "expression_interpretation_status": status,
                "phase7e_expression_support_level": support_level,
                "artifact_risk": artifact_risk,
                "conservative_interpretation": conservative_interpretation(
                    status,
                    gene,
                    len(detected),
                    len(quantified),
                    median(tpm_values) if tpm_values else None,
                    low_mapping_count,
                    very_low_mapping_count,
                ),
                "required_validation": required_validation,
                "claim_language_guardrail": "Use Phase 7e only for retina-specific detected/not-detected expression support. Do not infer activation, differential expression, pathway state, organism-wide aging, causation, functional advantage, gene absence, or longevity mechanism.",
                "supporting_files": combine_artifact_risk(
                    first["supporting_files"],
                    "results/rnaseq/phase7d_candidate_expression_matrix.tsv",
                    "results/rnaseq/phase7d_salmon_quant_summary.tsv",
                    "results/rnaseq/phase7c_candidate_reference_validation.tsv",
                    "results/orthology/candidate_duplication_audit.tsv",
                    phase4e.get("supporting_files", NOT_ASSESSED) if phase4e else NOT_ASSESSED,
                    phase5c.get("supporting_files", NOT_ASSESSED),
                ),
            }
        )
    return rows


def build_parameter_review_rows(config: dict, run_qc_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    phase7d = config["phase7d_rnaseq_quantification"]
    phase7e = config["phase7e_expression_hardening"]
    mapping_rates = [parse_float(row["salmon_percent_mapped"]) for row in run_qc_rows]
    mapping_rates = [value for value in mapping_rates if value is not None]
    decoys = join_values([row["salmon_decoy_count"] for row in run_qc_rows])
    library_types = join_values([row["salmon_library_types"] for row in run_qc_rows])
    rows = [
        {
            "parameter": "salmon_algorithm",
            "configured_value": "Salmon quant with --validateMappings",
            "observed_value": phase7d["salmon_executable"],
            "review_status": "SUITABLE_FOR_FIRST_PASS_CANDIDATE_QUANTIFICATION",
            "technical_rationale": "Salmon is efficient for transcript-level paired-end RNA-seq quantification and handles multi-mapping probabilistically, which is appropriate for a first-pass candidate matrix.",
            "scientific_guardrail": "Transcript-level quantification does not resolve paralog identity for duplicated or repeat-context loci.",
            "recommended_action": "Use genome-aligned counting or targeted locus-specific review for FTH1B, H1F0, RAD51, TP53, and other ambiguous candidates before strong Phase 8 weighting.",
        },
        {
            "parameter": "salmon_library_type",
            "configured_value": str(phase7d.get("salmon_library_type", NOT_ASSESSED)),
            "observed_value": library_types,
            "review_status": "SUITABLE_AUTO_DETECTION_FOR_CURRENT_FIRST_PASS",
            "technical_rationale": "The configured A mode lets Salmon infer orientation; all completed runs report ISR.",
            "scientific_guardrail": "Auto-detection should be recorded and reviewed, not hidden.",
            "recommended_action": "For a confirmatory rerun, consider pinning ISR after verifying library metadata and Salmon logs.",
        },
        {
            "parameter": "salmon_bias_flags",
            "configured_value": ";".join(str(arg) for arg in phase7d.get("salmon_extra_args", [])),
            "observed_value": "Quant logs completed with configured flags.",
            "review_status": "SUITABLE_FOR_FIRST_PASS",
            "technical_rationale": "--seqBias and --gcBias reduce common sequence/composition biases; --validateMappings improves mapping specificity compared with lightweight quasi-mapping alone.",
            "scientific_guardrail": "Bias correction and selective alignment do not make candidate expression causal or mechanism-level evidence.",
            "recommended_action": "Keep these flags for consistency; record any future changes as a new quantification run.",
        },
        {
            "parameter": "salmon_index_decoys",
            "configured_value": str(phase7d["salmon_index_dir"]),
            "observed_value": decoys,
            "review_status": "FIRST_PASS_ACCEPTABLE_CONFIRMATORY_DECOY_INDEX_RECOMMENDED" if decoys == "0" else "DECOY_STATUS_RECORDED",
            "technical_rationale": "The current index produced usable quant.sf files, but no-decoy indexing can be less conservative for off-transcriptome or paralogous sequence contexts.",
            "scientific_guardrail": "No-decoy Salmon results should be treated as first-pass candidate support only.",
            "recommended_action": "Build a decoy-aware index or run genome-aligned validation before using expression as strong Phase 8 support for ambiguous genes.",
        },
        {
            "parameter": "mapping_rate_thresholds",
            "configured_value": f"warning<{phase7e['low_mapping_rate_warning_percent']}%; very_low<{phase7e['very_low_mapping_rate_percent']}%",
            "observed_value": f"min={min(mapping_rates):.6f}; max={max(mapping_rates):.6f}" if mapping_rates else NOT_ASSESSED,
            "review_status": "SUITABLE_REVIEW_THRESHOLDS_NOT_STATISTICAL_FILTERS",
            "technical_rationale": "Thresholds flag runs for interpretation review without discarding completed quantifications.",
            "scientific_guardrail": "Low mapping rate is technical uncertainty, not biological evidence.",
            "recommended_action": "Carry mapping-rate caveats into Phase 8 and prioritize reference/index validation if expression support becomes important.",
        },
        {
            "parameter": "detection_thresholds",
            "configured_value": f"TPM>={phase7d['min_tpm_for_detected']}; NumReads>={phase7d['min_numreads_for_detected']}",
            "observed_value": "Applied in Phase 7d candidate matrix.",
            "review_status": "SUITABLE_EXPLORATORY_DETECTION_THRESHOLDS",
            "technical_rationale": "Requiring both TPM and read support avoids treating tiny TPM estimates as detected candidate expression.",
            "scientific_guardrail": "These thresholds are not a statistical expression model and cannot justify differential expression or pathway activity.",
            "recommended_action": "Keep for exploratory detected/not-detected summaries; use designed statistics only if future metadata support differential expression.",
        },
        {
            "parameter": "fastq_qc_sampling",
            "configured_value": str(phase7d["fastq_qc_max_records"]),
            "observed_value": "All six read files passed lightweight QC.",
            "review_status": "SUITABLE_SANITY_QC_NOT_PUBLICATION_GRADE_QC",
            "technical_rationale": "Sampling 100000 records per read is fast and catches malformed files, read length, base composition, and coarse quality issues.",
            "scientific_guardrail": "It does not replace full FastQC/MultiQC or contamination assessment.",
            "recommended_action": "Run full read QC if Phase 7 expression evidence becomes a central result.",
        },
        {
            "parameter": "minimum_consistent_retina_support",
            "configured_value": f"detected_runs>={phase7e['min_detected_runs_for_consistent_retina_support']}; fraction>={phase7e['min_detection_fraction_for_consistent_retina_support']}",
            "observed_value": f"quantified_runs={len(run_qc_rows)}",
            "review_status": "SUITABLE_FOR_THREE_RUN_RETINA_AUDIT",
            "technical_rationale": "For three retina runs, requiring at least two detected runs prevents one-run-only signals from being treated as consistent tissue support.",
            "scientific_guardrail": "This is not differential expression and does not generalize beyond retina.",
            "recommended_action": "Carry only consistent, low-ambiguity retina support forward as cautious Phase 8 expression evidence.",
        },
    ]
    return rows


def write_report(
    run_rows: list[dict[str, str]],
    candidate_rows: list[dict[str, str]],
    parameter_rows: list[dict[str, str]],
    output: Path,
) -> None:
    very_low = sum(1 for row in run_rows if row["mapping_rate_status"] == "VERY_LOW_MAPPING_RATE_REVIEW_REQUIRED")
    low = sum(1 for row in run_rows if row["mapping_rate_status"] == "LOW_MAPPING_RATE_CAUTIOUS_USE")
    phase8_ready = sum(1 for row in candidate_rows if row["phase7e_expression_support_level"] == "CAUTIOUS_RETINA_EXPRESSION_SUPPORT_FOR_PHASE8")
    deferred = sum(1 for row in candidate_rows if row["phase7e_expression_support_level"] == "DEFER_PHASE8_EXPRESSION_SUPPORT_UNTIL_LOCUS_REVIEW")
    lines = [
        "# Phase 7e Expression Interpretation Hardening",
        "",
        "Phase 7e reviews Phase 7d Salmon quantification outputs before expression support is allowed into Phase 8 evidence scoring.",
        "",
        "## Technical Status",
        "",
        f"- Quantified retina RNA-seq runs reviewed: {len(run_rows)}",
        f"- Low-mapping runs flagged: {low}",
        f"- Very-low-mapping runs flagged: {very_low}",
        f"- Candidate genes ready only for cautious retina expression support: {phase8_ready}",
        f"- Candidate genes deferred until locus/reference review: {deferred}",
        "",
        "## Method Rationale",
        "",
        "The implemented method is deterministic table integration: it joins Salmon run metrics, FASTQ sanity QC, Phase 7c reference ambiguity, Phase 3/4 duplication and locus-review outputs, and Phase 5 repeat-context artifact flags. This is preferred over a new statistical model because the current design contains three retina runs without a defensible contrast for differential expression.",
        "",
        "Salmon remains appropriate for first-pass candidate quantification because it is fast, reproducible, and records transcript-level estimates, but ambiguous paralogous loci require genome-aware validation before interpretation.",
        "",
        "## Guardrail",
        "",
        "Phase 7e may support cautious retina-specific detected/not-detected language. It does not support activation, differential expression, pathway state, organism-wide aging interpretation, causation, functional advantage, gene absence, or longevity mechanism.",
        "",
        "Supporting tables:",
        "",
        "- `results/rnaseq/phase7e_run_qc_review.tsv`",
        "- `results/rnaseq/phase7e_candidate_expression_hardened.tsv`",
        "- `results/rnaseq/phase7e_parameter_review.tsv`",
    ]
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase7e(
    config_path: Path,
    run_qc_output: Path,
    candidate_output: Path,
    parameter_review_output: Path,
    report_output: Path,
) -> None:
    config = read_yaml(config_path)
    if "phase7e_expression_hardening" not in config:
        raise ValueError("config.yaml is missing phase7e_expression_hardening")
    phase = config["phase7e_expression_hardening"]
    quant_rows = read_tsv(Path(phase["phase7d_salmon_quant_summary"]), PHASE7D_SALMON_QUANT_SUMMARY_COLUMNS)
    fastq_qc_rows = read_tsv(Path(phase["phase7d_fastq_qc"]), PHASE7D_FASTQ_QC_COLUMNS)
    expression_rows = read_tsv(Path(phase["phase7d_candidate_expression_matrix"]), PHASE7D_CANDIDATE_EXPRESSION_MATRIX_COLUMNS)
    phase7c_rows = read_tsv(Path(phase["phase7c_candidate_reference_validation"]), PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS)
    duplication_rows = read_tsv(Path(phase["candidate_duplication_audit"]), DUPLICATION_AUDIT_COLUMNS)
    phase4e_rows = read_tsv(Path(phase["phase4e_gene_hardened_summary"]), PHASE4E_GENE_HARDENING_COLUMNS)
    phase5c_rows = read_tsv(Path(phase["phase5c_gene_repeat_qc_summary"]), PHASE5C_GENE_QC_COLUMNS)

    run_rows = build_run_qc_rows(config, quant_rows, fastq_qc_rows)
    candidate_rows = build_candidate_hardened_rows(
        config,
        expression_rows,
        phase7c_rows,
        duplication_rows,
        phase4e_rows,
        phase5c_rows,
        run_rows,
    )
    parameter_rows = build_parameter_review_rows(config, run_rows)

    write_tsv(run_qc_output, run_rows, PHASE7E_RUN_QC_COLUMNS)
    write_tsv(candidate_output, candidate_rows, PHASE7E_CANDIDATE_HARDENED_COLUMNS)
    write_tsv(parameter_review_output, parameter_rows, PHASE7E_PARAMETER_REVIEW_COLUMNS)
    write_report(run_rows, candidate_rows, parameter_rows, report_output)
    LOGGER.info("Wrote Phase 7e outputs for %d runs and %d candidate genes", len(run_rows), len(candidate_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Harden Phase 7 retina RNA-seq expression interpretation before Phase 8 scoring.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-qc-output", type=Path, required=True)
    parser.add_argument("--candidate-output", type=Path, required=True)
    parser.add_argument("--parameter-review-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_phase7e(args.config, args.run_qc_output, args.candidate_output, args.parameter_review_output, args.report_output)


if __name__ == "__main__":
    main()
