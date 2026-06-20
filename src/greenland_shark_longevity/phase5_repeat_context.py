"""Phase 5 limited repeat/transposon context around hardened candidate loci."""

from __future__ import annotations

import argparse
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from .evidence import validate_guardrail_language
from .schemas import (
    LOCAL_FILE_INVENTORY_COLUMNS,
    PHASE4E_LOCUS_HARDENING_COLUMNS,
    PHASE5_GENE_REPEAT_CONTEXT_COLUMNS,
    PHASE5_LOCUS_REPEAT_CONTEXT_COLUMNS,
    PHASE5_REPEAT_FEATURE_COLUMNS,
    PHASE5_REPEAT_RESOURCE_COLUMNS,
)
from .utils import read_tsv, read_yaml, write_tsv

LOGGER = logging.getLogger(__name__)
NOT_ASSESSED = "NOT_ASSESSED"

REPEAT_FEATURE_TYPES = {
    "repeat_region",
    "mobile_element",
    "transposable_element",
    "transposable_element_gene",
    "ltr_retrotransposon",
    "long_terminal_repeat",
    "line",
    "sine",
    "dna_transposon",
    "terminal_inverted_repeat_element",
    "target_site_duplication",
    "repeat_unit",
    "direct_repeat",
    "dispersed_repeat",
}

REPEAT_FILE_KEYWORDS = (
    "repeat",
    "repeatmask",
    "repeatmodel",
    "rmsk",
    "rmout",
    "transpos",
    "transposable",
    "retrotrans",
    "te_annotation",
)

REPEAT_ATTRIBUTE_KEYS = (
    "repeat_class",
    "repeat_family",
    "repeat_name",
    "rpt_type",
    "rpt_family",
    "rpt_unit_seq",
)


@dataclass(frozen=True)
class RepeatFeature:
    repeat_id: str
    seqid: str
    start: int
    end: int
    strand: str
    repeat_class: str
    repeat_family: str
    repeat_name: str
    source_file: str
    parse_status: str
    notes: str

    def row(self) -> dict[str, str]:
        return {
            "repeat_id": self.repeat_id,
            "seqid": self.seqid,
            "start": str(self.start),
            "end": str(self.end),
            "strand": self.strand,
            "repeat_class": self.repeat_class,
            "repeat_family": self.repeat_family,
            "repeat_name": self.repeat_name,
            "source_file": self.source_file,
            "parse_status": self.parse_status,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class CandidateWindow:
    seqid: str
    start: int
    end: int


def parse_int(value: str | None, default: int = 0) -> int:
    if value in {None, "", NOT_ASSESSED}:
        return default
    return int(float(str(value)))


def parse_float(value: str | None, default: float = 0.0) -> float:
    if value in {None, "", NOT_ASSESSED}:
        return default
    return float(str(value))


def parse_gff_attributes(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in value.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, attr_value = part.split("=", 1)
            attrs[key.strip()] = unquote(attr_value.strip())
        elif " " in part:
            key, attr_value = part.split(" ", 1)
            attrs[key.strip()] = attr_value.strip().strip('"')
    return attrs


def interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    left = max(start_a, start_b)
    right = min(end_a, end_b)
    return max(0, right - left + 1)


def interval_distance(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    if interval_overlap(start_a, end_a, start_b, end_b):
        return 0
    if end_a < start_b:
        return start_b - end_a - 1
    return start_a - end_b - 1


def normalize_seqid(seqid: str, seqid_aliases: dict[str, str] | None = None) -> str:
    aliases = seqid_aliases or {}
    if seqid in aliases:
        return aliases[seqid]
    versionless = seqid.rsplit(".", 1)[0] if "." in seqid else seqid
    return aliases.get(versionless, seqid)


def load_assembly_seqid_aliases(assembly_report: Path | None) -> dict[str, str]:
    """Map NCBI assembly-report accessions to local sequence names."""
    if not assembly_report or not str(assembly_report) or not assembly_report.exists():
        return {}
    aliases: dict[str, str] = {}
    with assembly_report.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            sequence_name = fields[0]
            genbank_accn = fields[4]
            refseq_accn = fields[6]
            ucsc_name = fields[9]
            for value in (sequence_name, genbank_accn, refseq_accn, ucsc_name):
                if value and value != "na":
                    aliases[value] = sequence_name
                    if "." in value:
                        aliases[value.rsplit(".", 1)[0]] = sequence_name
            aliases[sequence_name] = sequence_name
    return aliases


def build_candidate_windows(
    candidate_loci_path: Path,
    window_bp: int,
    seqid_aliases: dict[str, str] | None = None,
) -> dict[str, list[CandidateWindow]]:
    windows: dict[str, list[CandidateWindow]] = {}
    if not candidate_loci_path.exists():
        return windows
    for locus in read_tsv(candidate_loci_path, PHASE4E_LOCUS_HARDENING_COLUMNS):
        start = parse_int(locus["locus_start"])
        end = parse_int(locus["locus_end"])
        window = CandidateWindow(
            seqid=normalize_seqid(locus["annotation_seqid"], seqid_aliases),
            start=max(1, start - window_bp),
            end=end + window_bp,
        )
        windows.setdefault(window.seqid, []).append(window)
        target_seqid = locus.get("target_seqid", "")
        if target_seqid:
            windows.setdefault(normalize_seqid(target_seqid, seqid_aliases), []).append(window)
    return windows


def feature_in_candidate_windows(
    seqid: str,
    start: int,
    end: int,
    target_windows: dict[str, list[CandidateWindow]] | None = None,
) -> bool:
    if not target_windows:
        return True
    return any(interval_overlap(start, end, window.start, window.end) for window in target_windows.get(seqid, []))


def has_repeat_filename_signal(path: Path) -> bool:
    lowered = path.name.lower()
    return any(keyword in lowered for keyword in REPEAT_FILE_KEYWORDS)


def detected_repeat_format(path: Path, inventory_format: str = "") -> str:
    suffixes = "".join(path.suffixes).lower()
    if inventory_format in {"gff3", "gff", "bed", "repeatmasker_out"}:
        return inventory_format
    if suffixes.endswith((".gff3", ".gff", ".gtf")):
        return "gff3"
    if suffixes.endswith(".bed"):
        return "bed"
    if path.name.lower().endswith((".out", ".rm.out")) and has_repeat_filename_signal(path):
        return "repeatmasker_out"
    return inventory_format or "UNKNOWN"


def gff_line_is_repeat(fields: list[str], attrs: dict[str, str], candidate_by_name: bool) -> bool:
    feature_type = fields[2].lower()
    if feature_type in REPEAT_FEATURE_TYPES:
        return True
    lowered_attrs = {key.lower(): value.lower() for key, value in attrs.items()}
    if any(key in lowered_attrs for key in REPEAT_ATTRIBUTE_KEYS):
        return True
    if candidate_by_name and fields[1].lower() == "repeatmasker":
        return True
    if candidate_by_name and lowered_attrs.get("target", "").startswith("motif:"):
        return True
    if candidate_by_name and any(keyword in " ".join(lowered_attrs.values()) for keyword in REPEAT_FILE_KEYWORDS):
        return True
    return False


def repeat_class_from_attrs(feature_type: str, attrs: dict[str, str]) -> str:
    for key in ("repeat_class", "class", "rpt_type", "mobile_element_type"):
        if attrs.get(key):
            return attrs[key]
    if feature_type.lower() in REPEAT_FEATURE_TYPES:
        return feature_type
    return NOT_ASSESSED


def repeat_family_from_attrs(attrs: dict[str, str]) -> str:
    for key in ("repeat_family", "family", "rpt_family"):
        if attrs.get(key):
            return attrs[key]
    repeat_class = attrs.get("class", "")
    if "/" in repeat_class:
        return repeat_class.split("/", 1)[1]
    return NOT_ASSESSED


def repeat_name_from_attrs(attrs: dict[str, str], fallback: str) -> str:
    for key in ("Name", "name", "repeat_name", "Target", "ID"):
        if attrs.get(key):
            value = attrs[key]
            target_token = value.split()[0].strip('"')
            if target_token.startswith("Motif:"):
                return target_token.replace("Motif:", "").strip('"')
            return value
    return fallback


def parse_gff_repeat_features(
    path: Path,
    candidate_by_name: bool = False,
    seqid_aliases: dict[str, str] | None = None,
    target_windows: dict[str, list[CandidateWindow]] | None = None,
) -> list[RepeatFeature]:
    features: list[RepeatFeature] = []
    generated = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = parse_gff_attributes(fields[8])
            if not gff_line_is_repeat(fields, attrs, candidate_by_name):
                continue
            original_seqid = fields[0]
            seqid = normalize_seqid(original_seqid, seqid_aliases)
            start = parse_int(fields[3])
            end = parse_int(fields[4])
            if not feature_in_candidate_windows(seqid, start, end, target_windows):
                continue
            generated += 1
            repeat_id = attrs.get("ID", f"{path.stem}_repeat_{generated:06d}")
            notes = "Parsed from GFF/GFF3 repeat-like feature or repeat attributes."
            if original_seqid != seqid:
                notes += f" Original seqid {original_seqid} normalized to {seqid}."
            if target_windows:
                notes += " Retained by candidate-window filter."
            features.append(
                RepeatFeature(
                    repeat_id=repeat_id,
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=fields[6],
                    repeat_class=repeat_class_from_attrs(fields[2], attrs),
                    repeat_family=repeat_family_from_attrs(attrs),
                    repeat_name=repeat_name_from_attrs(attrs, repeat_id),
                    source_file=str(path),
                    parse_status="PARSED_REPEAT_FEATURE",
                    notes=notes,
                )
            )
    return features


def parse_bed_repeat_features(
    path: Path,
    seqid_aliases: dict[str, str] | None = None,
    target_windows: dict[str, list[CandidateWindow]] | None = None,
) -> list[RepeatFeature]:
    features: list[RepeatFeature] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for idx, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            repeat_name = fields[3] if len(fields) > 3 and fields[3] else f"{path.stem}_bed_{idx:06d}"
            strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-", "."} else "."
            original_seqid = fields[0]
            seqid = normalize_seqid(original_seqid, seqid_aliases)
            start = parse_int(fields[1]) + 1
            end = parse_int(fields[2])
            if not feature_in_candidate_windows(seqid, start, end, target_windows):
                continue
            notes = "BED coordinates were converted from 0-based start to 1-based inclusive start."
            if original_seqid != seqid:
                notes += f" Original seqid {original_seqid} normalized to {seqid}."
            if target_windows:
                notes += " Retained by candidate-window filter."
            features.append(
                RepeatFeature(
                    repeat_id=f"{path.stem}_bed_{idx:06d}",
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    repeat_class=fields[6] if len(fields) > 6 and fields[6] else NOT_ASSESSED,
                    repeat_family=fields[7] if len(fields) > 7 and fields[7] else NOT_ASSESSED,
                    repeat_name=repeat_name,
                    source_file=str(path),
                    parse_status="PARSED_BED_INTERVAL",
                    notes=notes,
                )
            )
    return features


def parse_repeatmasker_out(
    path: Path,
    seqid_aliases: dict[str, str] | None = None,
    target_windows: dict[str, list[CandidateWindow]] | None = None,
) -> list[RepeatFeature]:
    features: list[RepeatFeature] = []
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
            strand = "-" if fields[8] == "C" else "+"
            repeat_name = fields[9]
            class_family = fields[10]
            repeat_class, repeat_family = (class_family.split("/", 1) + [NOT_ASSESSED])[:2] if "/" in class_family else (class_family, NOT_ASSESSED)
            notes = "Parsed from RepeatMasker .out coordinates."
            if original_seqid != seqid:
                notes += f" Original seqid {original_seqid} normalized to {seqid}."
            if target_windows:
                notes += " Retained by candidate-window filter."
            features.append(
                RepeatFeature(
                    repeat_id=f"{path.stem}_{len(features) + 1:06d}",
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    repeat_class=repeat_class,
                    repeat_family=repeat_family,
                    repeat_name=repeat_name,
                    source_file=str(path),
                    parse_status="PARSED_REPEATMASKER_OUT",
                    notes=notes,
                )
            )
    return features


def parse_repeat_features(
    path: Path,
    detected_format: str = "",
    candidate_by_name: bool = False,
    seqid_aliases: dict[str, str] | None = None,
    target_windows: dict[str, list[CandidateWindow]] | None = None,
) -> list[RepeatFeature]:
    fmt = detected_repeat_format(path, detected_format)
    if fmt == "gff3":
        return parse_gff_repeat_features(
            path,
            candidate_by_name=candidate_by_name,
            seqid_aliases=seqid_aliases,
            target_windows=target_windows,
        )
    if fmt == "bed":
        return parse_bed_repeat_features(path, seqid_aliases=seqid_aliases, target_windows=target_windows)
    if fmt == "repeatmasker_out":
        return parse_repeatmasker_out(path, seqid_aliases=seqid_aliases, target_windows=target_windows)
    return []


def inspect_repeat_resource(
    row: dict[str, str],
    seqid_aliases: dict[str, str] | None = None,
    target_windows: dict[str, list[CandidateWindow]] | None = None,
) -> dict[str, str]:
    path = Path(row["local_path"])
    detected_format = detected_repeat_format(path, row.get("detected_format", ""))
    base = {
        "resource_id": row.get("resource_id", "SMIC_TOKYO_GENOME_2025"),
        "source_id": row.get("source_id", "NOT_ASSESSED"),
        "file_name": row.get("file_name", path.name),
        "file_role_candidate": row.get("file_role_candidate", NOT_ASSESSED),
        "detected_format": detected_format,
        "local_path": str(path),
        "byte_size": row.get("byte_size", str(path.stat().st_size) if path.exists() else NOT_ASSESSED),
        "repeat_annotation_status": "NOT_REPEAT_ANNOTATION_CANDIDATE",
        "repeat_feature_count": NOT_ASSESSED,
        "inspection_method": "filename_and_format_screen",
        "notes": "File does not have a repeat-annotation filename or supported repeat-annotation format signal.",
    }
    if not path.exists():
        base.update(
            {
                "repeat_annotation_status": "MISSING_LOCAL_FILE",
                "notes": "Configured local path was not found.",
            }
        )
        return base
    if detected_format == "zip":
        base.update(
            {
                "repeat_annotation_status": "ARCHIVE_RETAINED_NOT_PARSED",
                "notes": "Source archive is retained for provenance; extracted contents are inspected separately.",
            }
        )
        return base
    if detected_format == "gff3":
        candidate_by_name = has_repeat_filename_signal(path)
        features = parse_gff_repeat_features(
            path,
            candidate_by_name=candidate_by_name,
            seqid_aliases=seqid_aliases,
            target_windows=target_windows,
        )
        filter_note = " Candidate-window filtering and sequence-id normalization were applied." if target_windows else ""
        base.update(
            {
                "repeat_annotation_status": "REPEAT_ANNOTATION_AVAILABLE" if features else "NO_REPEAT_FEATURES_DETECTED",
                "repeat_feature_count": str(len(features)),
                "inspection_method": "gff_repeat_feature_scan",
                "notes": (
                    f"GFF/GFF3 contains parseable repeat-like features.{filter_note}"
                    if features
                    else f"GFF/GFF3 was scanned for repeat-like feature types and repeat attributes; none were detected in the inspected scope.{filter_note}"
                ),
            }
        )
        return base
    if has_repeat_filename_signal(path) and detected_format in {"bed", "repeatmasker_out"}:
        features = parse_repeat_features(
            path,
            detected_format=detected_format,
            candidate_by_name=True,
            seqid_aliases=seqid_aliases,
            target_windows=target_windows,
        )
        filter_note = " Candidate-window filtering and sequence-id normalization were applied." if target_windows else ""
        base.update(
            {
                "repeat_annotation_status": "REPEAT_ANNOTATION_AVAILABLE" if features else "NO_PARSEABLE_REPEAT_FEATURES",
                "repeat_feature_count": str(len(features)),
                "inspection_method": f"{detected_format}_repeat_feature_scan",
                "notes": f"Repeat-annotation filename signal was found and the supported interval format was parsed.{filter_note}",
            }
        )
        return base
    return base


def explicit_repeat_rows(paths: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path_value in paths:
        path = Path(path_value)
        rows.append(
            {
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "source_id": "CONFIG_EXPLICIT_REPEAT_ANNOTATION",
                "file_name": path.name,
                "file_role_candidate": "explicit_repeat_annotation",
                "detected_format": detected_repeat_format(path),
                "local_path": str(path),
                "byte_size": str(path.stat().st_size) if path.exists() else NOT_ASSESSED,
                "sequence_count": "NOT_APPLICABLE",
                "notes": "Explicit repeat annotation candidate from config/config.yaml.",
            }
        )
    return rows


def discover_repeat_resources(
    figshare_inventory_path: Path,
    explicit_paths: list[str] | None = None,
    seqid_aliases: dict[str, str] | None = None,
    target_windows: dict[str, list[CandidateWindow]] | None = None,
) -> list[dict[str, str]]:
    candidate_rows: list[dict[str, str]] = []
    if figshare_inventory_path.exists():
        candidate_rows.extend(read_tsv(figshare_inventory_path, LOCAL_FILE_INVENTORY_COLUMNS))
    candidate_rows.extend(explicit_repeat_rows(explicit_paths or []))

    resource_rows = [inspect_repeat_resource(row, seqid_aliases=seqid_aliases, target_windows=target_windows) for row in candidate_rows]
    if not any(row["repeat_annotation_status"] == "REPEAT_ANNOTATION_AVAILABLE" for row in resource_rows):
        resource_rows.append(
            {
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "source_id": "PHASE5_LOCAL_REPEAT_SEARCH",
                "file_name": "NOT_APPLICABLE",
                "file_role_candidate": "local_repeat_annotation_search_summary",
                "detected_format": NOT_ASSESSED,
                "local_path": str(figshare_inventory_path),
                "byte_size": NOT_ASSESSED,
                "repeat_annotation_status": "NO_REPEAT_ANNOTATION_FOUND",
                "repeat_feature_count": "0",
                "inspection_method": "figshare_inventory_and_gff_repeat_feature_scan",
                "notes": "No parseable repeat annotation file was found in the current local Tokyo/PNAS Figshare inventory or explicit Phase 5 config paths after applying configured candidate-window filters.",
            }
        )
    return resource_rows


def load_repeat_features(
    resource_rows: list[dict[str, str]],
    seqid_aliases: dict[str, str] | None = None,
    target_windows: dict[str, list[CandidateWindow]] | None = None,
) -> list[RepeatFeature]:
    features: list[RepeatFeature] = []
    for row in resource_rows:
        if row["repeat_annotation_status"] != "REPEAT_ANNOTATION_AVAILABLE":
            continue
        path = Path(row["local_path"])
        features.extend(
            parse_repeat_features(
                path,
                row["detected_format"],
                candidate_by_name=has_repeat_filename_signal(path),
                seqid_aliases=seqid_aliases,
                target_windows=target_windows,
            )
        )
    return sorted(features, key=lambda item: (item.seqid, item.start, item.end, item.repeat_id))


def join_unique(values: list[str], empty: str = NOT_ASSESSED) -> str:
    cleaned = sorted({value for value in values if value and value != NOT_ASSESSED})
    return ";".join(cleaned) if cleaned else empty


def repeat_annotation_status(resource_rows: list[dict[str, str]]) -> str:
    if any(row["repeat_annotation_status"] == "REPEAT_ANNOTATION_AVAILABLE" for row in resource_rows):
        return "REPEAT_ANNOTATION_AVAILABLE"
    return "NO_REPEAT_ANNOTATION_AVAILABLE"


def build_locus_context_rows(
    candidate_loci_path: Path,
    repeat_features: list[RepeatFeature],
    resource_rows: list[dict[str, str]],
    window_bp: int,
) -> list[dict[str, str]]:
    loci = read_tsv(candidate_loci_path, PHASE4E_LOCUS_HARDENING_COLUMNS)
    status = repeat_annotation_status(resource_rows)
    by_seqid: dict[str, list[RepeatFeature]] = {}
    for feature in repeat_features:
        by_seqid.setdefault(feature.seqid, []).append(feature)

    rows: list[dict[str, str]] = []
    for locus in loci:
        start = parse_int(locus["locus_start"])
        end = parse_int(locus["locus_end"])
        window_start = max(1, start - window_bp)
        window_end = end + window_bp
        span = max(1, end - start + 1)
        seqid_features = by_seqid.get(locus["annotation_seqid"], [])

        if status != "REPEAT_ANNOTATION_AVAILABLE":
            row = {
                "gene_symbol": locus["gene_symbol"],
                "mechanism": locus["mechanism"],
                "locus_cluster_id": locus["locus_cluster_id"],
                "annotation_seqid": locus["annotation_seqid"],
                "locus_start": str(start),
                "locus_end": str(end),
                "locus_span_bp": str(span),
                "window_start": str(window_start),
                "window_end": str(window_end),
                "window_size_bp": str(window_end - window_start + 1),
                "repeat_annotation_status": status,
                "repeat_overlap_count": NOT_ASSESSED,
                "repeat_overlap_bp": NOT_ASSESSED,
                "repeat_overlap_fraction": NOT_ASSESSED,
                "nearest_repeat_distance_bp": NOT_ASSESSED,
                "repeat_classes": NOT_ASSESSED,
                "repeat_families": NOT_ASSESSED,
                "repeat_names": NOT_ASSESSED,
                "artifact_context_status": "REPEAT_CONTEXT_NOT_ASSESSABLE_CURRENT_LOCAL_FILES",
                "biological_interpretation": "Repeat-context evidence is not assessable because no local repeat annotation file was detected in the current Tokyo/PNAS/Figshare package.",
                "claim_language_guardrail": "Do not interpret Phase 5 local repeat context as biological mechanism, validated duplication, adaptation, pathway activity, or longevity evidence.",
                "required_validation": "Obtain or generate versioned repeat annotations, then rerun local interval intersection before interpreting repeat context.",
                "supporting_files": "data/metadata/phase5_repeat_resource_status.tsv;results/rescue/phase4e_locus_manual_review.tsv",
            }
            rows.append(row)
            continue

        overlapping = [feature for feature in seqid_features if interval_overlap(start, end, feature.start, feature.end)]
        window_features = [feature for feature in seqid_features if interval_overlap(window_start, window_end, feature.start, feature.end)]
        overlap_bp = sum(interval_overlap(start, end, feature.start, feature.end) for feature in overlapping)
        nearest = min((interval_distance(start, end, feature.start, feature.end) for feature in seqid_features), default=None)
        context_status = (
            "REPEAT_OVERLAP_RECORDED_ARTIFACT_CONTEXT_ONLY"
            if overlapping
            else ("LOCAL_REPEAT_FEATURES_IN_WINDOW_ONLY" if window_features else "NO_LOCAL_REPEAT_OVERLAP_IN_AVAILABLE_ANNOTATION")
        )
        rows.append(
            {
                "gene_symbol": locus["gene_symbol"],
                "mechanism": locus["mechanism"],
                "locus_cluster_id": locus["locus_cluster_id"],
                "annotation_seqid": locus["annotation_seqid"],
                "locus_start": str(start),
                "locus_end": str(end),
                "locus_span_bp": str(span),
                "window_start": str(window_start),
                "window_end": str(window_end),
                "window_size_bp": str(window_end - window_start + 1),
                "repeat_annotation_status": status,
                "repeat_overlap_count": str(len(overlapping)),
                "repeat_overlap_bp": str(overlap_bp),
                "repeat_overlap_fraction": f"{overlap_bp / span:.6g}",
                "nearest_repeat_distance_bp": str(nearest) if nearest is not None else NOT_ASSESSED,
                "repeat_classes": join_unique([feature.repeat_class for feature in window_features], empty="NONE_IN_WINDOW"),
                "repeat_families": join_unique([feature.repeat_family for feature in window_features], empty="NONE_IN_WINDOW"),
                "repeat_names": join_unique([feature.repeat_name for feature in window_features], empty="NONE_IN_WINDOW"),
                "artifact_context_status": context_status,
                "biological_interpretation": "Local repeat context was recorded as artifact/context evidence only and is not interpreted as functional or longevity evidence.",
                "claim_language_guardrail": "Do not interpret Phase 5 local repeat context as biological mechanism, validated duplication, adaptation, pathway activity, or longevity evidence.",
                "required_validation": "Confirm repeat annotation provenance, repeat class assignments, candidate-locus coordinates, and cross-resource support before using repeat context in Phase 8 scoring.",
                "supporting_files": "results/repeats/phase5_repeat_features.tsv;results/rescue/phase4e_locus_manual_review.tsv",
            }
        )
    validate_guardrail_language(
        [
            {
                "gene_or_pathway": row["gene_symbol"],
                "biological_interpretation": row["biological_interpretation"],
                "relevance_to_aging_longevity": NOT_ASSESSED,
                "translational_relevance": NOT_ASSESSED,
            }
            for row in rows
        ]
    )
    return rows


def build_gene_context_rows(locus_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_gene: dict[str, list[dict[str, str]]] = {}
    for row in locus_rows:
        by_gene.setdefault(row["gene_symbol"], []).append(row)

    rows: list[dict[str, str]] = []
    for gene_symbol, gene_rows in sorted(by_gene.items()):
        status = gene_rows[0]["repeat_annotation_status"]
        starts = [parse_int(row["locus_start"]) for row in gene_rows]
        ends = [parse_int(row["locus_end"]) for row in gene_rows]
        seqids = join_unique([row["annotation_seqid"] for row in gene_rows])
        if status != "REPEAT_ANNOTATION_AVAILABLE":
            rows.append(
                {
                    "gene_symbol": gene_symbol,
                    "mechanism": gene_rows[0]["mechanism"],
                    "locus_count": str(len(gene_rows)),
                    "annotation_seqids": seqids,
                    "cluster_start": str(min(starts)),
                    "cluster_end": str(max(ends)),
                    "repeat_annotation_status": status,
                    "loci_with_repeat_overlap": NOT_ASSESSED,
                    "total_repeat_overlap_bp": NOT_ASSESSED,
                    "max_locus_repeat_overlap_fraction": NOT_ASSESSED,
                    "repeat_classes": NOT_ASSESSED,
                    "repeat_families": NOT_ASSESSED,
                    "phase5_context_status": "REPEAT_CONTEXT_NOT_ASSESSABLE_CURRENT_LOCAL_FILES",
                    "artifact_risk": NOT_ASSESSED,
                    "conservative_interpretation": "Repeat context is not assessable from the current local Tokyo/PNAS/Figshare files because no repeat annotation file was detected.",
                    "required_validation": "Use versioned repeat annotations or a documented optional repeat-discovery workflow before making repeat-context interpretations.",
                    "claim_language_guardrail": "Do not interpret Phase 5 local repeat context as biological mechanism, validated duplication, adaptation, pathway activity, or longevity evidence.",
                    "supporting_files": "data/metadata/phase5_repeat_resource_status.tsv;results/rescue/phase4e_locus_manual_review.tsv",
                }
            )
            continue

        overlap_counts = [parse_int(row["repeat_overlap_count"]) for row in gene_rows]
        overlap_bp = [parse_int(row["repeat_overlap_bp"]) for row in gene_rows]
        overlap_fractions = [parse_float(row["repeat_overlap_fraction"]) for row in gene_rows]
        rows.append(
            {
                "gene_symbol": gene_symbol,
                "mechanism": gene_rows[0]["mechanism"],
                "locus_count": str(len(gene_rows)),
                "annotation_seqids": seqids,
                "cluster_start": str(min(starts)),
                "cluster_end": str(max(ends)),
                "repeat_annotation_status": status,
                "loci_with_repeat_overlap": str(sum(1 for value in overlap_counts if value > 0)),
                "total_repeat_overlap_bp": str(sum(overlap_bp)),
                "max_locus_repeat_overlap_fraction": f"{max(overlap_fractions) if overlap_fractions else 0.0:.6g}",
                "repeat_classes": join_unique([row["repeat_classes"] for row in gene_rows if row["repeat_classes"] != "NONE_IN_WINDOW"], empty="NONE_IN_WINDOW"),
                "repeat_families": join_unique([row["repeat_families"] for row in gene_rows if row["repeat_families"] != "NONE_IN_WINDOW"], empty="NONE_IN_WINDOW"),
                "phase5_context_status": "LOCAL_REPEAT_CONTEXT_RECORDED_ARTIFACT_CONTEXT_ONLY",
                "artifact_risk": "context_recorded_not_interpreted",
                "conservative_interpretation": "Local repeat annotations were intersected with candidate loci as artifact/context evidence only; this does not support functional or longevity interpretation.",
                "required_validation": "Confirm repeat annotation provenance, repeat classifications, locus coordinates, and cross-resource support before Phase 8 scoring.",
                "claim_language_guardrail": "Do not interpret Phase 5 local repeat context as biological mechanism, validated duplication, adaptation, pathway activity, or longevity evidence.",
                "supporting_files": "results/repeats/phase5_candidate_locus_repeat_context.tsv;results/repeats/phase5_repeat_features.tsv",
            }
        )
    validate_guardrail_language(
        [
            {
                "gene_or_pathway": row["gene_symbol"],
                "biological_interpretation": row["conservative_interpretation"],
                "relevance_to_aging_longevity": NOT_ASSESSED,
                "translational_relevance": NOT_ASSESSED,
            }
            for row in rows
        ]
    )
    return rows


def write_phase5_outputs(
    config_path: Path,
    candidate_loci: Path,
    figshare_inventory: Path,
    resource_status_output: Path,
    repeat_features_output: Path,
    locus_context_output: Path,
    gene_summary_output: Path,
) -> None:
    config = read_yaml(config_path) if config_path.exists() else {}
    phase_config = config.get("phase5_repeat_context", {})
    explicit_paths = [str(path) for path in phase_config.get("repeat_annotation_candidates", [])]
    window_bp = parse_int(str(phase_config.get("candidate_window_bp", 10000)), default=10000)
    assembly_report = Path(str(phase_config.get("assembly_report", ""))) if phase_config.get("assembly_report") else None
    filter_to_candidate_windows = bool(phase_config.get("filter_to_candidate_windows", False))
    seqid_aliases = load_assembly_seqid_aliases(assembly_report)
    target_windows = build_candidate_windows(candidate_loci, window_bp, seqid_aliases) if filter_to_candidate_windows else None

    if seqid_aliases:
        LOGGER.info("Loaded %d sequence-id aliases from %s", len(seqid_aliases), assembly_report)
    if target_windows:
        LOGGER.info("Applying candidate-window repeat filtering across %d normalized seqids", len(target_windows))

    resource_rows = discover_repeat_resources(
        figshare_inventory,
        explicit_paths,
        seqid_aliases=seqid_aliases,
        target_windows=target_windows,
    )
    features = load_repeat_features(resource_rows, seqid_aliases=seqid_aliases, target_windows=target_windows)
    locus_rows = build_locus_context_rows(candidate_loci, features, resource_rows, window_bp)
    gene_rows = build_gene_context_rows(locus_rows)

    write_tsv(resource_status_output, resource_rows, PHASE5_REPEAT_RESOURCE_COLUMNS)
    write_tsv(repeat_features_output, [feature.row() for feature in features], PHASE5_REPEAT_FEATURE_COLUMNS)
    write_tsv(locus_context_output, locus_rows, PHASE5_LOCUS_REPEAT_CONTEXT_COLUMNS)
    write_tsv(gene_summary_output, gene_rows, PHASE5_GENE_REPEAT_CONTEXT_COLUMNS)
    LOGGER.info("Wrote Phase 5 repeat context for %d repeat features and %d candidate loci", len(features), len(locus_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 5 limited repeat context for candidate loci.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--candidate-loci", type=Path, required=True)
    parser.add_argument("--figshare-inventory", type=Path, required=True)
    parser.add_argument("--resource-status-output", type=Path, required=True)
    parser.add_argument("--repeat-features-output", type=Path, required=True)
    parser.add_argument("--locus-context-output", type=Path, required=True)
    parser.add_argument("--gene-summary-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    write_phase5_outputs(
        args.config,
        args.candidate_loci,
        args.figshare_inventory,
        args.resource_status_output,
        args.repeat_features_output,
        args.locus_context_output,
        args.gene_summary_output,
    )


if __name__ == "__main__":
    main()
