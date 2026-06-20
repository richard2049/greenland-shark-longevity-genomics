"""Phase 7d retina RNA-seq raw-read intake, QC, and Salmon quantification."""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
from pathlib import Path

from .schemas import (
    PHASE7_RNASEQ_MANIFEST_COLUMNS,
    PHASE7_RNASEQ_READINESS_COLUMNS,
    PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS,
    PHASE7C_TX2GENE_COLUMNS,
    PHASE7D_CANDIDATE_EXPRESSION_MATRIX_COLUMNS,
    PHASE7D_FASTQ_QC_COLUMNS,
    PHASE7D_RAW_READ_INTAKE_COLUMNS,
    PHASE7D_SALMON_PREFLIGHT_COLUMNS,
    PHASE7D_SALMON_QUANT_SUMMARY_COLUMNS,
)
from .utils import NOT_ASSESSED, ensure_parent, join_values, open_text, read_tsv, read_yaml, write_tsv

LOGGER = logging.getLogger(__name__)


def bool_from_config(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def format_float(value: float | None, digits: int = 4) -> str:
    if value is None:
        return NOT_ASSESSED
    return f"{value:.{digits}f}"


def index_by_run(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["run"]: row for row in rows}


def eligible_retina_runs(manifest_rows: list[dict[str, str]], readiness_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    readiness_by_run = index_by_run(readiness_rows)
    rows: list[dict[str, str]] = []
    for row in manifest_rows:
        readiness = readiness_by_run.get(row["run"])
        if readiness is None:
            continue
        if readiness.get("included_for_candidate_expression_plan") != "true":
            continue
        if row.get("library_strategy") != "RNA-Seq":
            continue
        if row.get("tissue", "").lower() != "retina":
            continue
        rows.append(row)
    return rows


def expected_fastq_paths(run: str, raw_dir: Path, patterns: list[str]) -> tuple[Path, Path]:
    candidates = [raw_dir / pattern.format(run=run) for pattern in patterns]
    read1_candidates = [path for path in candidates if "_1." in path.name or "_R1." in path.name]
    read2_candidates = [path for path in candidates if "_2." in path.name or "_R2." in path.name]
    read1 = next((path for path in read1_candidates if path.exists()), read1_candidates[0])
    read2 = next((path for path in read2_candidates if path.exists()), read2_candidates[0])
    return read1, read2


def file_status(path: Path) -> str:
    if path.exists() and path.is_file() and path.stat().st_size > 0:
        return "LOCAL_AVAILABLE"
    if path.exists() and path.is_file():
        return "EMPTY_FILE"
    return "MISSING_LOCAL"


def build_raw_read_intake_rows(config: dict, manifest_rows: list[dict[str, str]], readiness_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    phase = config["phase7d_rnaseq_quantification"]
    raw_dir = Path(phase["raw_read_local_dir"])
    patterns = phase["expected_fastq_patterns"]
    rows: list[dict[str, str]] = []
    for run_row in eligible_retina_runs(manifest_rows, readiness_rows):
        read1, read2 = expected_fastq_paths(run_row["run"], raw_dir, patterns)
        read1_status = file_status(read1)
        read2_status = file_status(read2)
        pair_ready = read1_status == "LOCAL_AVAILABLE" and read2_status == "LOCAL_AVAILABLE"
        rows.append(
            {
                "resource_id": run_row["resource_id"],
                "run": run_row["run"],
                "experiment": run_row["experiment"],
                "biosample": run_row["biosample"],
                "sample_name": run_row["sample_name"],
                "tissue": run_row["tissue"],
                "library_layout": run_row["library_layout"],
                "read1_path": str(read1),
                "read2_path": str(read2),
                "read1_status": read1_status,
                "read2_status": read2_status,
                "read1_size_bytes": str(read1.stat().st_size) if read1.exists() else "0",
                "read2_size_bytes": str(read2.stat().st_size) if read2.exists() else "0",
                "raw_pair_status": "PAIRED_FASTQ_READY" if pair_ready else "PAIRED_FASTQ_MISSING_OR_INCOMPLETE",
                "md5_status": "NOT_ASSESSED_CONFIG_CHECKSUM_RAW_READS_FALSE",
                "intake_status": "READY_FOR_FASTQ_QC" if pair_ready else "RAW_READS_NOT_AVAILABLE_LOCALLY",
                "required_action": "Run FASTQ QC and Salmon quantification." if pair_ready else "Download or link paired FASTQ files for this retina RNA-seq run, then rerun Phase 7d.",
            }
        )
    if not rows:
        raise ValueError("No retina RNA-seq runs are eligible for Phase 7d")
    return rows


def inspect_fastq(path: Path, max_records: int) -> dict[str, str]:
    if file_status(path) != "LOCAL_AVAILABLE":
        return {
            "file_status": file_status(path),
            "records_inspected": "0",
            "mean_read_length": NOT_ASSESSED,
            "gc_percent": NOT_ASSESSED,
            "n_percent": NOT_ASSESSED,
            "mean_phred_quality": NOT_ASSESSED,
            "qc_status": "NOT_RUN_FASTQ_MISSING",
            "notes": "FASTQ file is not available locally.",
        }
    records = 0
    total_bases = 0
    gc = 0
    n_count = 0
    quality_sum = 0
    quality_bases = 0
    malformed = False
    with open_text(path) as handle:
        while records < max_records:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip().upper()
            plus = handle.readline()
            quality = handle.readline().strip()
            if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                malformed = True
                break
            records += 1
            total_bases += len(sequence)
            gc += sequence.count("G") + sequence.count("C")
            n_count += sequence.count("N")
            quality_sum += sum(ord(char) - 33 for char in quality)
            quality_bases += len(quality)
    if malformed:
        return {
            "file_status": "LOCAL_AVAILABLE",
            "records_inspected": str(records),
            "mean_read_length": NOT_ASSESSED,
            "gc_percent": NOT_ASSESSED,
            "n_percent": NOT_ASSESSED,
            "mean_phred_quality": NOT_ASSESSED,
            "qc_status": "FAIL_MALFORMED_FASTQ_RECORD",
            "notes": "A malformed FASTQ record was encountered during lightweight QC.",
        }
    if records == 0 or total_bases == 0:
        return {
            "file_status": "LOCAL_AVAILABLE",
            "records_inspected": str(records),
            "mean_read_length": NOT_ASSESSED,
            "gc_percent": NOT_ASSESSED,
            "n_percent": NOT_ASSESSED,
            "mean_phred_quality": NOT_ASSESSED,
            "qc_status": "FAIL_NO_READS_INSPECTED",
            "notes": "No complete FASTQ records were inspected.",
        }
    return {
        "file_status": "LOCAL_AVAILABLE",
        "records_inspected": str(records),
        "mean_read_length": format_float(total_bases / records, 2),
        "gc_percent": format_float((gc / total_bases) * 100, 3),
        "n_percent": format_float((n_count / total_bases) * 100, 3),
        "mean_phred_quality": format_float(quality_sum / quality_bases, 3) if quality_bases else NOT_ASSESSED,
        "qc_status": "PASS_LIGHTWEIGHT_FASTQ_QC",
        "notes": f"QC inspected up to {max_records} reads; this is not a substitute for full FastQC/MultiQC when publication-grade read QC is required.",
    }


def build_fastq_qc_rows(intake_rows: list[dict[str, str]], max_records: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in intake_rows:
        for read_label, path_column in [("R1", "read1_path"), ("R2", "read2_path")]:
            path = Path(row[path_column])
            qc = inspect_fastq(path, max_records)
            rows.append(
                {
                    "run": row["run"],
                    "read": read_label,
                    "fastq_path": str(path),
                    **qc,
                }
            )
    return rows


def salmon_path(executable: str) -> str | None:
    path = Path(executable)
    if path.exists():
        return str(path)
    return shutil.which(executable)


def index_ready(index_dir: Path) -> bool:
    return index_dir.exists() and any(index_dir.iterdir())


def build_salmon_preflight_rows(config: dict, intake_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    phase = config["phase7d_rnaseq_quantification"]
    salmon = salmon_path(phase["salmon_executable"])
    transcript_fasta = Path(phase["transcript_fasta"])
    tx2gene = Path(phase["phase7c_tx2gene"])
    index_dir = Path(phase["salmon_index_dir"])
    all_fastqs_ready = all(row["raw_pair_status"] == "PAIRED_FASTQ_READY" for row in intake_rows)
    rows = [
        {
            "check_id": "salmon_executable",
            "item": "Salmon executable",
            "path_or_value": salmon or phase["salmon_executable"],
            "status": "PASS" if salmon else "MISSING",
            "message": "Salmon executable is available." if salmon else "Salmon is not available on PATH or configured executable path.",
            "required_action": "None" if salmon else "Install Salmon in WSL/Linux, Docker, or the active environment before quantification.",
        },
        {
            "check_id": "transcript_fasta",
            "item": "Phase 7c transcript FASTA",
            "path_or_value": str(transcript_fasta),
            "status": "PASS" if transcript_fasta.exists() else "MISSING",
            "message": "Transcript FASTA is available." if transcript_fasta.exists() else "Phase 7c transcript FASTA is missing.",
            "required_action": "None" if transcript_fasta.exists() else "Run Phase 7c before Salmon indexing/quantification.",
        },
        {
            "check_id": "tx2gene",
            "item": "Phase 7c tx2gene",
            "path_or_value": str(tx2gene),
            "status": "PASS" if tx2gene.exists() else "MISSING",
            "message": "tx2gene table is available." if tx2gene.exists() else "Phase 7c tx2gene table is missing.",
            "required_action": "None" if tx2gene.exists() else "Run Phase 7c before candidate summarization.",
        },
        {
            "check_id": "salmon_index",
            "item": "Salmon index directory",
            "path_or_value": str(index_dir),
            "status": "PASS" if index_ready(index_dir) else "MISSING",
            "message": "Salmon index is available." if index_ready(index_dir) else "Salmon index is missing or empty.",
            "required_action": "None" if index_ready(index_dir) else "Build a Salmon index from the Phase 7c transcript FASTA with recorded parameters.",
        },
        {
            "check_id": "paired_fastqs",
            "item": "Eligible retina paired FASTQ files",
            "path_or_value": join_values([row["run"] for row in intake_rows if row["raw_pair_status"] == "PAIRED_FASTQ_READY"]),
            "status": "PASS" if all_fastqs_ready else "MISSING",
            "message": "All eligible retina RNA-seq FASTQ pairs are available." if all_fastqs_ready else "One or more eligible retina RNA-seq FASTQ pairs are missing.",
            "required_action": "None" if all_fastqs_ready else "Download or link paired FASTQ files, then rerun Phase 7d.",
        },
    ]
    return rows


def build_index_if_requested(config: dict, salmon: str | None, preflight_rows: list[dict[str, str]]) -> None:
    phase = config["phase7d_rnaseq_quantification"]
    if not bool_from_config(phase.get("build_salmon_index_if_missing", False)):
        return
    index_dir = Path(phase["salmon_index_dir"])
    transcript_fasta = Path(phase["transcript_fasta"])
    if salmon is None or index_ready(index_dir) or not transcript_fasta.exists():
        return
    ensure_parent(index_dir / "placeholder")
    command = [salmon, "index", "-t", str(transcript_fasta), "-i", str(index_dir)]
    LOGGER.info("Running Salmon index: %s", " ".join(command))
    subprocess.run(command, check=True, text=True)
    for row in preflight_rows:
        if row["check_id"] == "salmon_index":
            row["status"] = "PASS" if index_ready(index_dir) else "MISSING"
            row["message"] = "Salmon index was built by Phase 7d." if index_ready(index_dir) else row["message"]


def salmon_quant_command(config: dict, run: str, read1: Path, read2: Path, output_dir: Path) -> list[str]:
    phase = config["phase7d_rnaseq_quantification"]
    salmon = salmon_path(phase["salmon_executable"])
    if salmon is None:
        raise ValueError("Salmon executable is not available")
    command = [
        salmon,
        "quant",
        "-i",
        str(Path(phase["salmon_index_dir"])),
        "-l",
        str(phase.get("salmon_library_type", "A")),
        "-1",
        str(read1),
        "-2",
        str(read2),
        "-p",
        str(phase.get("salmon_threads", 1)),
        "-o",
        str(output_dir),
    ]
    command.extend(str(arg) for arg in phase.get("salmon_extra_args", []))
    return command


def run_salmon_quantification(config: dict, intake_rows: list[dict[str, str]], preflight_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    phase = config["phase7d_rnaseq_quantification"]
    salmon = salmon_path(phase["salmon_executable"])
    run_salmon = bool_from_config(phase.get("run_salmon_if_ready", True))
    reuse_existing = bool_from_config(phase.get("reuse_existing_salmon_quant_if_present", True))
    index_dir = Path(phase["salmon_index_dir"])
    quant_root = Path(phase["salmon_quant_dir"])
    logs_dir = Path(phase["logs_dir"])
    rows: list[dict[str, str]] = []
    for row in intake_rows:
        output_dir = quant_root / row["run"]
        quant_sf = output_dir / "quant.sf"
        stdout_log = logs_dir / f"{row['run']}.salmon.stdout.log"
        stderr_log = logs_dir / f"{row['run']}.salmon.stderr.log"
        read1 = Path(row["read1_path"])
        read2 = Path(row["read2_path"])
        command_text = "NOT_RUN"
        status = "NOT_RUN"
        message = ""
        if reuse_existing and quant_sf.exists() and quant_sf.stat().st_size > 0:
            command_text = "REUSED_EXISTING_QUANT_SF"
            status = "SALMON_QUANT_COMPLETE"
            message = "Existing quant.sf reused as a restart checkpoint; delete the sample Salmon output directory to force re-quantification."
        elif row["raw_pair_status"] != "PAIRED_FASTQ_READY":
            status = "NOT_RUN_MISSING_RAW_FASTQ_PAIR"
            message = "Paired FASTQ files are not available locally."
        elif salmon is None:
            status = "NOT_RUN_MISSING_SALMON_EXECUTABLE"
            message = "Salmon executable is unavailable."
        elif not index_ready(index_dir):
            status = "NOT_RUN_MISSING_SALMON_INDEX"
            message = "Salmon index is unavailable."
        elif not run_salmon:
            status = "NOT_RUN_CONFIG_RUN_SALMON_FALSE"
            message = "Config disables Salmon execution."
        else:
            ensure_parent(stdout_log)
            ensure_parent(output_dir / "placeholder")
            command = salmon_quant_command(config, row["run"], read1, read2, output_dir)
            command_text = " ".join(command)
            LOGGER.info("Running Salmon quant for %s", row["run"])
            with stdout_log.open("w", encoding="utf-8") as stdout, stderr_log.open("w", encoding="utf-8") as stderr:
                subprocess.run(command, check=True, stdout=stdout, stderr=stderr, text=True)
            status = "SALMON_QUANT_COMPLETE" if quant_sf.exists() else "SALMON_QUANT_FINISHED_MISSING_QUANT_SF"
            message = "Salmon quantification completed." if quant_sf.exists() else "Salmon command finished but quant.sf was not found."
        rows.append(
            {
                "run": row["run"],
                "tissue": row["tissue"],
                "read1_path": str(read1),
                "read2_path": str(read2),
                "salmon_output_dir": str(output_dir),
                "quant_sf": str(quant_sf),
                "quantification_status": status,
                "salmon_command": command_text,
                "log_stdout": str(stdout_log),
                "log_stderr": str(stderr_log),
                "message": message,
            }
        )
    return rows


def parse_quant_sf(path: Path) -> dict[str, dict[str, float]]:
    if not path.exists():
        return {}
    rows = read_tsv(path, ["Name", "TPM", "NumReads"])
    parsed: dict[str, dict[str, float]] = {}
    for row in rows:
        parsed[row["Name"]] = {
            "TPM": float(row["TPM"]),
            "NumReads": float(row["NumReads"]),
        }
    return parsed


def expression_status(tpm: float, numreads: float, min_tpm: float, min_numreads: float, quant_status: str, readiness: str) -> str:
    if quant_status != "SALMON_QUANT_COMPLETE":
        return "NOT_QUANTIFIED"
    if readiness.startswith("REFERENCE_PRESENT_TARGETED") or readiness.startswith("REFERENCE_PRESENT_BUT_NOT_READY"):
        return "QUANTIFIED_REQUIRES_LOCUS_REVIEW_BEFORE_INTERPRETATION"
    if tpm >= min_tpm and numreads >= min_numreads:
        return "DETECTED_IN_RETINA_EXPLORATORY"
    return "NOT_DETECTED_UNDER_CURRENT_RETINA_QUANTIFICATION"


def build_candidate_expression_matrix(
    config: dict,
    quant_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    phase = config["phase7d_rnaseq_quantification"]
    candidate_rows = read_tsv(Path(phase["phase7c_candidate_reference_validation"]), PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS)
    min_tpm = float(phase.get("min_tpm_for_detected", 1.0))
    min_numreads = float(phase.get("min_numreads_for_detected", 10.0))
    rows: list[dict[str, str]] = []
    for quant_row in quant_rows:
        quant = parse_quant_sf(Path(quant_row["quant_sf"]))
        for candidate in candidate_rows:
            transcript_ids = [value for value in candidate["matched_transcript_ids"].split(";") if value and value != NOT_ASSESSED]
            values = [quant[transcript_id] for transcript_id in transcript_ids if transcript_id in quant]
            summed_tpm = sum(value["TPM"] for value in values)
            summed_numreads = sum(value["NumReads"] for value in values)
            status = expression_status(
                summed_tpm,
                summed_numreads,
                min_tpm,
                min_numreads,
                quant_row["quantification_status"],
                candidate["phase7c_quantification_readiness"],
            )
            rows.append(
                {
                    "mechanism": candidate["mechanism"],
                    "gene_symbol": candidate["gene_symbol"],
                    "run": quant_row["run"],
                    "tissue": quant_row["tissue"],
                    "phase7c_quantification_readiness": candidate["phase7c_quantification_readiness"],
                    "matched_transcript_count": candidate["matched_transcript_count"],
                    "quantified_transcript_count": str(len(values)),
                    "summed_tpm": format_float(summed_tpm, 6) if quant_row["quantification_status"] == "SALMON_QUANT_COMPLETE" else NOT_ASSESSED,
                    "summed_numreads": format_float(summed_numreads, 3) if quant_row["quantification_status"] == "SALMON_QUANT_COMPLETE" else NOT_ASSESSED,
                    "candidate_expression_status": status,
                    "quantification_status": quant_row["quantification_status"],
                    "conservative_interpretation": "Retina-only candidate quantification summary; not differential expression and not pathway activity." if status != "NOT_QUANTIFIED" else "No expression interpretation because quantification was not run for this sample.",
                    "claim_language_guardrail": "Use detected/not detected in retina only after successful quantification. This is not differential expression. Do not infer activation, pathway state, causation, functional advantage, organism-wide aging, or longevity mechanism.",
                    "supporting_files": f"{quant_row['quant_sf']};{phase['phase7c_candidate_reference_validation']};{phase['phase7c_tx2gene']}",
                }
            )
    return rows


def write_report(
    intake_rows: list[dict[str, str]],
    fastq_qc_rows: list[dict[str, str]],
    quant_rows: list[dict[str, str]],
    matrix_rows: list[dict[str, str]],
    output: Path,
) -> None:
    ready_fastqs = sum(1 for row in intake_rows if row["raw_pair_status"] == "PAIRED_FASTQ_READY")
    qc_pass = sum(1 for row in fastq_qc_rows if row["qc_status"] == "PASS_LIGHTWEIGHT_FASTQ_QC")
    quant_complete = sum(1 for row in quant_rows if row["quantification_status"] == "SALMON_QUANT_COMPLETE")
    detected = sum(1 for row in matrix_rows if row["candidate_expression_status"] == "DETECTED_IN_RETINA_EXPLORATORY")
    lines = [
        "# Phase 7d Retina RNA-seq Quantification",
        "",
        "Phase 7d validates local paired FASTQ intake, runs lightweight FASTQ QC when files are present, and runs Salmon quantification only when all configured inputs are available.",
        "",
        "## Status",
        "",
        f"- Eligible retina RNA-seq runs: {len(intake_rows)}",
        f"- Runs with local paired FASTQs: {ready_fastqs}",
        f"- FASTQ files passing lightweight QC: {qc_pass}",
        f"- Salmon quantifications complete: {quant_complete}",
        f"- Candidate/run detected summaries: {detected}",
        "",
        "Guardrail: Phase 7d is retina-only. It does not perform differential expression and must not be interpreted as pathway activity, organism-wide aging state, causation, or longevity mechanism evidence.",
        "",
        "Supporting tables:",
        "",
        "- `results/rnaseq/phase7d_raw_read_intake.tsv`",
        "- `results/rnaseq/phase7d_fastq_qc.tsv`",
        "- `results/rnaseq/phase7d_salmon_preflight.tsv`",
        "- `results/rnaseq/phase7d_salmon_quant_summary.tsv`",
        "- `results/rnaseq/phase7d_candidate_expression_matrix.tsv`",
    ]
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase7d(
    config_path: Path,
    raw_read_intake_output: Path,
    fastq_qc_output: Path,
    salmon_preflight_output: Path,
    salmon_quant_summary_output: Path,
    candidate_expression_matrix_output: Path,
    report_output: Path,
) -> None:
    config = read_yaml(config_path)
    if "phase7d_rnaseq_quantification" not in config:
        raise ValueError("config.yaml is missing phase7d_rnaseq_quantification")
    phase = config["phase7d_rnaseq_quantification"]
    manifest_rows = read_tsv(Path(phase["phase7a_manifest"]), PHASE7_RNASEQ_MANIFEST_COLUMNS)
    readiness_rows = read_tsv(Path(phase["phase7a_readiness"]), PHASE7_RNASEQ_READINESS_COLUMNS)
    read_tsv(Path(phase["phase7c_tx2gene"]), PHASE7C_TX2GENE_COLUMNS)
    read_tsv(Path(phase["phase7c_candidate_reference_validation"]), PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS)

    intake_rows = build_raw_read_intake_rows(config, manifest_rows, readiness_rows)
    fastq_qc_rows = build_fastq_qc_rows(intake_rows, int(phase.get("fastq_qc_max_records", 100000)))
    preflight_rows = build_salmon_preflight_rows(config, intake_rows)
    build_index_if_requested(config, salmon_path(phase["salmon_executable"]), preflight_rows)
    quant_rows = run_salmon_quantification(config, intake_rows, preflight_rows)
    matrix_rows = build_candidate_expression_matrix(config, quant_rows)

    write_tsv(raw_read_intake_output, intake_rows, PHASE7D_RAW_READ_INTAKE_COLUMNS)
    write_tsv(fastq_qc_output, fastq_qc_rows, PHASE7D_FASTQ_QC_COLUMNS)
    write_tsv(salmon_preflight_output, preflight_rows, PHASE7D_SALMON_PREFLIGHT_COLUMNS)
    write_tsv(salmon_quant_summary_output, quant_rows, PHASE7D_SALMON_QUANT_SUMMARY_COLUMNS)
    write_tsv(candidate_expression_matrix_output, matrix_rows, PHASE7D_CANDIDATE_EXPRESSION_MATRIX_COLUMNS)
    write_report(intake_rows, fastq_qc_rows, quant_rows, matrix_rows, report_output)
    LOGGER.info("Wrote Phase 7d outputs for %d retina RNA-seq runs", len(intake_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Phase 7d retina RNA-seq intake/QC and guarded Salmon quantification.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--raw-read-intake-output", type=Path, required=True)
    parser.add_argument("--fastq-qc-output", type=Path, required=True)
    parser.add_argument("--salmon-preflight-output", type=Path, required=True)
    parser.add_argument("--salmon-quant-summary-output", type=Path, required=True)
    parser.add_argument("--candidate-expression-matrix-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_phase7d(
        args.config,
        args.raw_read_intake_output,
        args.fastq_qc_output,
        args.salmon_preflight_output,
        args.salmon_quant_summary_output,
        args.candidate_expression_matrix_output,
        args.report_output,
    )


if __name__ == "__main__":
    main()
