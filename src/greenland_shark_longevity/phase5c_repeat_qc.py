"""Phase 5c repeat-context QC hardening for candidate loci."""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path

from .schemas import (
    PHASE4E_LOCUS_HARDENING_COLUMNS,
    PHASE5C_GENE_QC_COLUMNS,
    PHASE5C_INTEGRITY_COLUMNS,
    PHASE5C_LOCUS_QC_COLUMNS,
    PHASE5B_OUTPUT_INVENTORY_COLUMNS,
    PHASE5_LOCUS_REPEAT_CONTEXT_COLUMNS,
)
from .phase5_repeat_context import (
    NOT_ASSESSED,
    CandidateWindow,
    build_candidate_windows,
    feature_in_candidate_windows,
    interval_overlap,
    load_assembly_seqid_aliases,
    normalize_seqid,
    parse_float,
    parse_int,
)
from .utils import ensure_parent, read_tsv, read_yaml, write_tsv

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class RepeatMaskerOutFeature:
    seqid: str
    start: int
    end: int
    strand: str
    repeat_name: str
    repeat_class: str
    repeat_family: str
    score: str
    perc_div: str


def split_class_family(value: str) -> tuple[str, str]:
    if "/" in value:
        repeat_class, repeat_family = value.split("/", 1)
        return repeat_class, repeat_family
    return value, NOT_ASSESSED


def join_unique(values: list[str], empty: str = NOT_ASSESSED) -> str:
    cleaned = sorted({value for value in values if value and value != NOT_ASSESSED})
    return ";".join(cleaned) if cleaned else empty


def clipped_interval(start: int, end: int, clip_start: int, clip_end: int) -> tuple[int, int] | None:
    clipped_start = max(start, clip_start)
    clipped_end = min(end, clip_end)
    if clipped_start > clipped_end:
        return None
    return clipped_start, clipped_end


def union_bp(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return sum(end - start + 1 for start, end in merged)


def parse_repeatmasker_out_features(
    path: Path,
    seqid_aliases: dict[str, str],
    target_windows: dict[str, list[CandidateWindow]],
) -> list[RepeatMaskerOutFeature]:
    """Stream RepeatMasker .out and retain only records intersecting candidate windows."""
    if not path.exists():
        return []
    features: list[RepeatMaskerOutFeature] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(("SW", "score", "There were")):
                continue
            fields = stripped.split()
            if len(fields) < 11 or not fields[0].lstrip("-").isdigit():
                continue
            original_seqid = fields[4]
            seqid = normalize_seqid(original_seqid, seqid_aliases)
            start = parse_int(fields[5])
            end = parse_int(fields[6])
            start, end = min(start, end), max(start, end)
            if not feature_in_candidate_windows(seqid, start, end, target_windows):
                continue
            repeat_class, repeat_family = split_class_family(fields[10])
            features.append(
                RepeatMaskerOutFeature(
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand="-" if fields[8] == "C" else "+",
                    repeat_name=fields[9],
                    repeat_class=repeat_class,
                    repeat_family=repeat_family,
                    score=fields[0],
                    perc_div=fields[1],
                )
            )
    LOGGER.info("Retained %d RepeatMasker .out records in candidate windows", len(features))
    return features


def load_phase5b_inventory(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    rows = read_tsv(path, PHASE5B_OUTPUT_INVENTORY_COLUMNS)
    inventory: dict[str, dict[str, str]] = {}
    for row in rows:
        inventory[row["output_id"]] = row
        inventory[row["local_path"]] = row
    return inventory


def log_warning_status(log_paths: list[Path]) -> tuple[str, str]:
    found_existing = False
    warnings: list[str] = []
    for path in log_paths:
        if not path.exists():
            continue
        found_existing = True
        text = path.read_text(encoding="utf-8", errors="replace")
        lowered = text.lower()
        if "signal 9" in lowered or "forksys" in lowered:
            warnings.append(f"{path}: ProcessRepeats signal-9/forksys warning detected")
        if "ltrpipeline" in lowered and ("no results" in lowered or "error" in lowered):
            warnings.append(f"{path}: LTRPipeline warning detected")
    if warnings:
        return "WARNING_DETECTED", "; ".join(warnings)
    if found_existing:
        return "NO_WARNING_DETECTED_IN_CHECKED_LOGS", "Checked configured logs for signal-9 and LTRPipeline warnings."
    return "LOG_NOT_AVAILABLE", "No configured RepeatMasker/RepeatModeler log was found."


def integrity_rows(config: dict, phase5b_inventory_path: Path) -> list[dict[str, str]]:
    phase5b = config.get("phase5b_repeat_annotation", {})
    phase5c = config.get("phase5c_repeat_qc", {})
    inventory = load_phase5b_inventory(phase5b_inventory_path)
    log_paths = [Path(path) for path in phase5c.get("log_files", ["logs/phase5b_repeatmasker.stdout.log"])]
    warning_status, warning_notes = log_warning_status(log_paths)

    roles = [
        ("repeatmodeler_library", phase5b.get("custom_library", "results/repeats/phase5b/smic_tokyo-families.fa")),
        ("repeatmasker_gff", phase5b.get("repeatmasker_gff", "results/repeats/phase5b/smic_tokyo.repeatmasker.out.gff")),
        ("repeatmasker_out", phase5b.get("repeatmasker_out", "results/repeats/phase5b/smic_tokyo.repeatmasker.out")),
        ("repeatmasker_tbl", phase5b.get("repeatmasker_tbl", "results/repeats/phase5b/smic_tokyo.repeatmasker.tbl")),
    ]
    rows: list[dict[str, str]] = []
    gff_present = Path(str(roles[1][1])).exists()
    out_present = Path(str(roles[2][1])).exists()
    tbl_present = Path(str(roles[3][1])).exists()
    if gff_present and out_present and warning_status == "WARNING_DETECTED":
        usability = "CANDIDATE_LOCUS_CONTEXT_ONLY_WITH_PROCESSREPEATS_WARNING"
    elif gff_present and out_present and tbl_present:
        usability = "REPEATMASKER_OUTPUT_AVAILABLE_FOR_CONTEXT_AND_SUMMARY"
    elif gff_present and out_present:
        usability = "CANDIDATE_LOCUS_CONTEXT_ONLY_NO_TBL_SUMMARY"
    else:
        usability = "REPEAT_CONTEXT_QC_NOT_ASSESSABLE"

    for idx, (role, path_value) in enumerate(roles, start=1):
        path = Path(str(path_value))
        inv = inventory.get(role) or inventory.get(str(path)) or {}
        status = "PRESENT" if path.exists() else "MISSING"
        notes = warning_notes
        if role == "repeatmasker_tbl" and not path.exists():
            notes += " RepeatMasker .tbl is missing, so genome-wide repeat summaries are not supported."
        rows.append(
            {
                "check_id": f"PHASE5C_INTEGRITY_{idx:03d}",
                "file_role": role,
                "path": str(path),
                "status": status,
                "byte_size": str(path.stat().st_size) if path.exists() else NOT_ASSESSED,
                "phase5b_inventory_status": inv.get("status", NOT_ASSESSED),
                "phase5b_inventory_sha256": inv.get("sha256", NOT_ASSESSED),
                "log_warning_status": warning_status,
                "candidate_context_usability": usability,
                "notes": notes,
            }
        )
    return rows


def context_by_locus(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path, PHASE5_LOCUS_REPEAT_CONTEXT_COLUMNS)
    return {row["locus_cluster_id"]: row for row in rows}


def out_features_by_seqid(features: list[RepeatMaskerOutFeature]) -> dict[str, list[RepeatMaskerOutFeature]]:
    by_seqid: dict[str, list[RepeatMaskerOutFeature]] = {}
    for feature in features:
        by_seqid.setdefault(feature.seqid, []).append(feature)
    return by_seqid


def artifact_modifier(gff_fraction: float, out_fraction: float, window_density: float) -> str:
    if max(gff_fraction, out_fraction) >= 0.5:
        return "high_repeat_context_artifact_risk"
    if max(gff_fraction, out_fraction) > 0:
        return "moderate_repeat_context_artifact_risk"
    if window_density > 0:
        return "repeat_context_window_only_no_direct_overlap"
    return "no_local_repeat_context_in_candidate_window"


def concordance_label(gff_count: int, out_count: int) -> str:
    if gff_count > 0 and out_count > 0:
        return "BOTH_GFF_AND_OUT_DIRECT_OVERLAP"
    if gff_count > 0:
        return "GFF_ONLY_DIRECT_OVERLAP_REVIEW_REQUIRED"
    if out_count > 0:
        return "OUT_ONLY_DIRECT_OVERLAP_REVIEW_REQUIRED"
    return "BOTH_GFF_AND_OUT_NO_DIRECT_OVERLAP"


def qc_status_from_concordance(concordance: str, window_density: float) -> str:
    if concordance == "BOTH_GFF_AND_OUT_DIRECT_OVERLAP":
        return "REPEAT_CONTEXT_SUPPORTED_BY_GFF_AND_OUT_ARTIFACT_ONLY"
    if "ONLY_DIRECT_OVERLAP" in concordance:
        return "REPEAT_CONTEXT_PARSER_DISCORDANCE_REVIEW_REQUIRED"
    if window_density > 0:
        return "WINDOW_REPEAT_CONTEXT_ONLY_NO_DIRECT_OVERLAP"
    return "NO_REPEAT_CONTEXT_IN_WINDOW_UNDER_CURRENT_REPEATMASKER_OUT"


def build_locus_qc_rows(
    candidate_loci_path: Path,
    phase5_locus_context_path: Path,
    out_features: list[RepeatMaskerOutFeature],
) -> list[dict[str, str]]:
    context_rows = context_by_locus(phase5_locus_context_path)
    by_seqid = out_features_by_seqid(out_features)
    rows: list[dict[str, str]] = []
    for locus in read_tsv(candidate_loci_path, PHASE4E_LOCUS_HARDENING_COLUMNS):
        context = context_rows.get(locus["locus_cluster_id"], {})
        seqid = locus["annotation_seqid"]
        start = parse_int(locus["locus_start"])
        end = parse_int(locus["locus_end"])
        window_start = parse_int(context.get("window_start"), max(1, start - 10000))
        window_end = parse_int(context.get("window_end"), end + 10000)
        window_size = max(1, window_end - window_start + 1)
        locus_span = max(1, end - start + 1)
        seqid_features = by_seqid.get(seqid, [])
        direct = [feature for feature in seqid_features if interval_overlap(start, end, feature.start, feature.end)]
        window = [feature for feature in seqid_features if interval_overlap(window_start, window_end, feature.start, feature.end)]
        direct_intervals = [
            interval
            for feature in direct
            for interval in [clipped_interval(feature.start, feature.end, start, end)]
            if interval is not None
        ]
        window_intervals = [
            interval
            for feature in window
            for interval in [clipped_interval(feature.start, feature.end, window_start, window_end)]
            if interval is not None
        ]
        out_overlap_bp = union_bp(direct_intervals)
        out_window_bp = union_bp(window_intervals)
        out_fraction = out_overlap_bp / locus_span
        window_density = out_window_bp / window_size
        gff_count = parse_int(context.get("repeat_overlap_count"))
        gff_fraction = parse_float(context.get("repeat_overlap_fraction"))
        concordance = concordance_label(gff_count, len(direct))
        rows.append(
            {
                "gene_symbol": locus["gene_symbol"],
                "mechanism": locus["mechanism"],
                "locus_cluster_id": locus["locus_cluster_id"],
                "annotation_seqid": seqid,
                "locus_start": str(start),
                "locus_end": str(end),
                "window_start": str(window_start),
                "window_end": str(window_end),
                "window_size_bp": str(window_size),
                "gff_repeat_overlap_count": str(gff_count),
                "gff_repeat_overlap_fraction": context.get("repeat_overlap_fraction", NOT_ASSESSED),
                "repeatmasker_out_overlap_count": str(len(direct)),
                "repeatmasker_out_overlap_bp": str(out_overlap_bp),
                "repeatmasker_out_overlap_fraction": f"{out_fraction:.6g}",
                "repeatmasker_out_window_repeat_bp": str(out_window_bp),
                "repeatmasker_out_window_repeat_density": f"{window_density:.6g}",
                "repeatmasker_out_classes": join_unique([feature.repeat_class for feature in window], empty="NONE_IN_WINDOW"),
                "repeatmasker_out_families": join_unique([feature.repeat_family for feature in window], empty="NONE_IN_WINDOW"),
                "repeatmasker_out_names": join_unique([feature.repeat_name for feature in window], empty="NONE_IN_WINDOW"),
                "gff_out_overlap_concordance": concordance,
                "qc_context_status": qc_status_from_concordance(concordance, window_density),
                "artifact_risk_modifier": artifact_modifier(gff_fraction, out_fraction, window_density),
                "interpretation_guardrail": "Use Phase 5c repeat context only as artifact risk; do not infer repeat-mediated mechanism, validated duplication, adaptation, pathway activity, causation, or longevity relevance.",
                "required_validation": "Manually review RepeatMasker GFF/.out intervals, candidate-locus coordinates, repeat class assignments, paralog identity, and cross-resource support before Phase 8 scoring.",
                "supporting_files": "results/repeats/phase5_candidate_locus_repeat_context.tsv;results/repeats/phase5b/smic_tokyo.repeatmasker.out;results/repeats/phase5c_locus_repeat_qc.tsv",
            }
        )
    return rows


def summarize_gene_qc(locus_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_gene: dict[str, list[dict[str, str]]] = {}
    for row in locus_rows:
        by_gene.setdefault(row["gene_symbol"], []).append(row)
    summaries: list[dict[str, str]] = []
    for gene_symbol, rows in sorted(by_gene.items()):
        out_overlap_fractions = [parse_float(row["repeatmasker_out_overlap_fraction"]) for row in rows]
        densities = [parse_float(row["repeatmasker_out_window_repeat_density"]) for row in rows]
        concordance_values = [row["gff_out_overlap_concordance"] for row in rows]
        artifact_values = [row["artifact_risk_modifier"] for row in rows]
        if "REPEAT_CONTEXT_PARSER_DISCORDANCE_REVIEW_REQUIRED" in {row["qc_context_status"] for row in rows}:
            qc_status = "REPEAT_CONTEXT_QC_DISCORDANT_REVIEW_REQUIRED"
        elif any(row["gff_out_overlap_concordance"] == "BOTH_GFF_AND_OUT_DIRECT_OVERLAP" for row in rows):
            qc_status = "REPEAT_CONTEXT_QC_SUPPORTED_BY_GFF_AND_OUT_ARTIFACT_ONLY"
        elif any(parse_float(row["repeatmasker_out_window_repeat_density"]) > 0 for row in rows):
            qc_status = "REPEAT_WINDOW_CONTEXT_QC_RECORDED_ARTIFACT_ONLY"
        else:
            qc_status = "NO_LOCAL_REPEAT_CONTEXT_IN_CURRENT_QC_WINDOW"
        if "high_repeat_context_artifact_risk" in artifact_values:
            artifact = "high_repeat_context_artifact_risk"
        elif "moderate_repeat_context_artifact_risk" in artifact_values:
            artifact = "moderate_repeat_context_artifact_risk"
        elif "repeat_context_window_only_no_direct_overlap" in artifact_values:
            artifact = "repeat_context_window_only_no_direct_overlap"
        else:
            artifact = "no_local_repeat_context_in_candidate_window"
        summaries.append(
            {
                "gene_symbol": gene_symbol,
                "mechanism": rows[0]["mechanism"],
                "locus_count": str(len(rows)),
                "annotation_seqids": join_unique([row["annotation_seqid"] for row in rows]),
                "gff_loci_with_direct_overlap": str(sum(1 for row in rows if parse_int(row["gff_repeat_overlap_count"]) > 0)),
                "repeatmasker_out_loci_with_direct_overlap": str(sum(1 for row in rows if parse_int(row["repeatmasker_out_overlap_count"]) > 0)),
                "max_repeatmasker_out_overlap_fraction": f"{max(out_overlap_fractions) if out_overlap_fractions else 0.0:.6g}",
                "max_window_repeat_density": f"{max(densities) if densities else 0.0:.6g}",
                "repeatmasker_out_classes": join_unique([row["repeatmasker_out_classes"] for row in rows if row["repeatmasker_out_classes"] != "NONE_IN_WINDOW"], empty="NONE_IN_WINDOW"),
                "repeatmasker_out_families": join_unique([row["repeatmasker_out_families"] for row in rows if row["repeatmasker_out_families"] != "NONE_IN_WINDOW"], empty="NONE_IN_WINDOW"),
                "gff_out_concordance_summary": join_unique(concordance_values),
                "phase5c_qc_status": qc_status,
                "artifact_risk_modifier": artifact,
                "conservative_interpretation": "Phase 5c compares GFF-derived repeat context with RepeatMasker .out intervals as artifact/context QC only; it does not support repeat-mediated biological interpretation.",
                "required_validation": "Review RepeatMasker interval provenance, local coordinates, repeat classes, paralog identity, and cross-resource consistency before Phase 8/Phase 9 interpretation.",
                "claim_language_guardrail": "Do not interpret Phase 5c repeat context as biological mechanism, validated duplication, adaptation, pathway activity, causation, translational relevance, or longevity evidence.",
                "supporting_files": "results/repeats/phase5c_locus_repeat_qc.tsv;results/repeats/phase5_candidate_locus_repeat_context.tsv;results/repeats/phase5b/smic_tokyo.repeatmasker.out",
            }
        )
    return summaries


def write_report(
    integrity: list[dict[str, str]],
    gene_rows: list[dict[str, str]],
    output: Path,
) -> None:
    usability = integrity[0]["candidate_context_usability"] if integrity else NOT_ASSESSED
    warning = integrity[0]["log_warning_status"] if integrity else NOT_ASSESSED
    lines = [
        "# Phase 5c Repeat-Context QC Hardening Report",
        "",
        f"RepeatMasker usability: `{usability}`.",
        f"Log warning status: `{warning}`.",
        "",
        "Phase 5c compares the GFF-derived candidate-locus repeat context against a window-filtered parse of the RepeatMasker `.out` file. This is an artifact-risk QC layer only.",
        "",
        "| Gene | QC status | GFF direct-overlap loci | `.out` direct-overlap loci | Max `.out` overlap fraction | Max window density | Artifact risk modifier |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for row in gene_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{row['gene_symbol']}`",
                    row["phase5c_qc_status"],
                    row["gff_loci_with_direct_overlap"],
                    row["repeatmasker_out_loci_with_direct_overlap"],
                    row["max_repeatmasker_out_overlap_fraction"],
                    row["max_window_repeat_density"],
                    row["artifact_risk_modifier"],
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "Interpretation guardrail: repeat overlap or density is treated as local artifact/context risk, not evidence of repeat-driven duplication, function, adaptation, or longevity.",
            "",
            "Supporting tables:",
            "",
            "- `results/repeats/phase5c_repeatmasker_integrity.tsv`",
            "- `results/repeats/phase5c_locus_repeat_qc.tsv`",
            "- `results/repeats/phase5c_gene_repeat_qc_summary.tsv`",
        ]
    )
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_phase5c_outputs(
    config_path: Path,
    candidate_loci: Path,
    phase5_locus_context: Path,
    phase5b_inventory: Path,
    integrity_output: Path,
    locus_qc_output: Path,
    gene_qc_output: Path,
    report_output: Path,
) -> None:
    config = read_yaml(config_path)
    phase5 = config.get("phase5_repeat_context", {})
    phase5b = config.get("phase5b_repeat_annotation", {})
    window_bp = parse_int(str(phase5.get("candidate_window_bp", 10000)), default=10000)
    assembly_report = Path(str(phase5.get("assembly_report", ""))) if phase5.get("assembly_report") else None
    seqid_aliases = load_assembly_seqid_aliases(assembly_report)
    target_windows = build_candidate_windows(candidate_loci, window_bp, seqid_aliases)
    repeatmasker_out = Path(str(phase5b.get("repeatmasker_out", "results/repeats/phase5b/smic_tokyo.repeatmasker.out")))

    integrity = integrity_rows(config, phase5b_inventory)
    out_features = parse_repeatmasker_out_features(repeatmasker_out, seqid_aliases, target_windows)
    locus_rows = build_locus_qc_rows(candidate_loci, phase5_locus_context, out_features)
    gene_rows = summarize_gene_qc(locus_rows)

    write_tsv(integrity_output, integrity, PHASE5C_INTEGRITY_COLUMNS)
    write_tsv(locus_qc_output, locus_rows, PHASE5C_LOCUS_QC_COLUMNS)
    write_tsv(gene_qc_output, gene_rows, PHASE5C_GENE_QC_COLUMNS)
    write_report(integrity, gene_rows, report_output)
    LOGGER.info("Wrote Phase 5c QC for %d candidate loci and %d genes", len(locus_rows), len(gene_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 5c repeat-context QC hardening.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--candidate-loci", type=Path, required=True)
    parser.add_argument("--phase5-locus-context", type=Path, required=True)
    parser.add_argument("--phase5b-inventory", type=Path, required=True)
    parser.add_argument("--integrity-output", type=Path, required=True)
    parser.add_argument("--locus-qc-output", type=Path, required=True)
    parser.add_argument("--gene-qc-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    write_phase5c_outputs(
        args.config,
        args.candidate_loci,
        args.phase5_locus_context,
        args.phase5b_inventory,
        args.integrity_output,
        args.locus_qc_output,
        args.gene_qc_output,
        args.report_output,
    )


if __name__ == "__main__":
    main()
