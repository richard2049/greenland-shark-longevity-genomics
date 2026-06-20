"""Phase 7c expression-reference construction and input validation."""

from __future__ import annotations

import argparse
import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from urllib.parse import unquote

from .schemas import (
    PHASE7B_CANDIDATE_QUANT_MAP_COLUMNS,
    PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS,
    PHASE7C_REFERENCE_QC_COLUMNS,
    PHASE7C_TX2GENE_COLUMNS,
)
from .utils import NOT_ASSESSED, ensure_parent, join_values, open_text, read_tsv, read_yaml, split_delimited_values, write_tsv

LOGGER = logging.getLogger(__name__)
SCAFFOLD_RE = re.compile(r"\bscaffold_\d+\b")
COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    gene_symbol: str
    locus_tag: str
    scaffold: str
    strand: str
    transcript_feature_type: str
    product: str = NOT_ASSESSED
    intervals: list[tuple[int, int]] = field(default_factory=list)
    interval_feature_type: str = "exon"


def parse_attributes(attribute_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for raw_part in attribute_text.split(";"):
        part = raw_part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def split_ids(value: str | None) -> list[str]:
    if not value:
        return []
    return [part.strip() for part in value.split(",") if part.strip()]


def scaffold_aliases(header: str) -> set[str]:
    text = header[1:].strip() if header.startswith(">") else header.strip()
    first_token = text.split()[0] if text else ""
    aliases = {first_token} if first_token else set()
    aliases.update(SCAFFOLD_RE.findall(text))
    return aliases


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1].upper()


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def parse_gff_models(
    annotation_gff: Path,
    transcript_feature_types: set[str],
    interval_feature: str,
    fallback_interval_feature: str,
) -> dict[str, TranscriptModel]:
    gene_metadata: dict[str, dict[str, str]] = {}
    transcripts: dict[str, TranscriptModel] = {}
    preferred_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    fallback_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)

    with annotation_gff.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if raw_line.startswith("#") or not raw_line.strip():
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) != 9:
                raise ValueError(f"Malformed GFF row in {annotation_gff} line {line_number}")
            scaffold, _source, feature_type, start, end, _score, strand, _phase, attributes = parts
            attrs = parse_attributes(attributes)
            if feature_type == "gene":
                gene_id = attrs.get("ID", NOT_ASSESSED)
                if gene_id != NOT_ASSESSED:
                    gene_metadata[gene_id] = attrs
                continue
            if feature_type in transcript_feature_types:
                transcript_id = attrs.get("ID", NOT_ASSESSED)
                if transcript_id == NOT_ASSESSED:
                    continue
                gene_id = attrs.get("Parent", NOT_ASSESSED)
                transcripts[transcript_id] = TranscriptModel(
                    transcript_id=transcript_id,
                    gene_id=gene_id,
                    gene_symbol=attrs.get("gene") or attrs.get("Name") or NOT_ASSESSED,
                    locus_tag=attrs.get("locus_tag", NOT_ASSESSED),
                    scaffold=scaffold,
                    strand=strand,
                    transcript_feature_type=feature_type,
                    product=attrs.get("product", NOT_ASSESSED),
                )
                continue
            if feature_type in {interval_feature, fallback_interval_feature}:
                target = preferred_intervals if feature_type == interval_feature else fallback_intervals
                for parent in split_ids(attrs.get("Parent")):
                    target[parent].append((int(start), int(end)))

    for transcript_id, model in transcripts.items():
        if model.gene_symbol == NOT_ASSESSED and model.gene_id in gene_metadata:
            gene_attrs = gene_metadata[model.gene_id]
            model.gene_symbol = gene_attrs.get("gene") or gene_attrs.get("Name") or NOT_ASSESSED
        if model.locus_tag == NOT_ASSESSED and model.gene_id in gene_metadata:
            model.locus_tag = gene_metadata[model.gene_id].get("locus_tag", NOT_ASSESSED)
        if transcript_id in preferred_intervals:
            model.intervals = sorted(preferred_intervals[transcript_id])
            model.interval_feature_type = interval_feature
        elif transcript_id in fallback_intervals:
            model.intervals = sorted(fallback_intervals[transcript_id])
            model.interval_feature_type = f"{fallback_interval_feature}_FALLBACK"

    usable = {key: model for key, model in transcripts.items() if model.intervals}
    if not usable:
        raise ValueError(f"No transcript models with {interval_feature} or {fallback_interval_feature} intervals found in {annotation_gff}")
    return usable


def format_fasta_header(model: TranscriptModel, length: int) -> str:
    gene_symbol = model.gene_symbol.replace(" ", "_")
    product = model.product.replace(" ", "_").replace(";", ",")
    return (
        f">{model.transcript_id} gene_id={model.gene_id} gene={gene_symbol} "
        f"locus_tag={model.locus_tag} scaffold={model.scaffold} strand={model.strand} "
        f"feature={model.interval_feature_type} exon_count={len(model.intervals)} "
        f"length={length} product={product}"
    )


def extract_transcript_sequence(scaffold_sequence: str, model: TranscriptModel) -> str:
    pieces = [scaffold_sequence[start - 1 : end] for start, end in sorted(model.intervals)]
    sequence = "".join(pieces).upper()
    if model.strand == "-":
        return reverse_complement(sequence)
    return sequence


def write_transcript_reference(
    genome_fasta: Path,
    transcript_fasta_output: Path,
    tx2gene_output: Path,
    transcript_models: dict[str, TranscriptModel],
) -> tuple[list[dict[str, str]], dict[str, int]]:
    models_by_scaffold: dict[str, list[TranscriptModel]] = defaultdict(list)
    for model in transcript_models.values():
        models_by_scaffold[model.scaffold].append(model)

    tx2gene_rows: list[dict[str, str]] = []
    stats = {
        "missing_scaffold_count": 0,
        "zero_length_transcript_count": 0,
        "exon_derived_transcript_count": 0,
        "cds_fallback_transcript_count": 0,
        "total_transcript_bp": 0,
    }
    seen_scaffolds: set[str] = set()

    ensure_parent(transcript_fasta_output)
    ensure_parent(tx2gene_output)
    with open_text(genome_fasta) as genome_handle, transcript_fasta_output.open("w", encoding="utf-8") as fasta_handle:
        current_header: str | None = None
        current_parts: list[str] = []

        def flush_record() -> None:
            if current_header is None:
                return
            aliases = scaffold_aliases(current_header)
            matched_scaffold = next((alias for alias in aliases if alias in models_by_scaffold), None)
            if matched_scaffold is None:
                return
            seen_scaffolds.add(matched_scaffold)
            scaffold_sequence = "".join(current_parts).upper()
            for model in models_by_scaffold[matched_scaffold]:
                sequence = extract_transcript_sequence(scaffold_sequence, model)
                if not sequence:
                    stats["zero_length_transcript_count"] += 1
                    continue
                length = len(sequence)
                stats["total_transcript_bp"] += length
                if model.interval_feature_type == "exon":
                    stats["exon_derived_transcript_count"] += 1
                else:
                    stats["cds_fallback_transcript_count"] += 1
                fasta_handle.write(format_fasta_header(model, length) + "\n")
                fasta_handle.write(wrap_sequence(sequence) + "\n")
                tx2gene_rows.append(
                    {
                        "transcript_id": model.transcript_id,
                        "gene_id": model.gene_id,
                        "gene_symbol": model.gene_symbol,
                        "locus_tag": model.locus_tag,
                        "scaffold": model.scaffold,
                        "strand": model.strand,
                        "feature_type": model.interval_feature_type,
                        "exon_count": str(len(model.intervals)),
                        "transcript_length_bp": str(length),
                        "reference_status": "TRANSCRIPT_REFERENCE_CONSTRUCTED",
                    }
                )

        for raw_line in genome_handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_record()
                current_header = line
                current_parts = []
            else:
                current_parts.append(line)
        flush_record()

    missing = set(models_by_scaffold) - seen_scaffolds
    stats["missing_scaffold_count"] = len(missing)
    if not tx2gene_rows:
        raise ValueError(f"No transcript sequences were written from {genome_fasta}")
    write_tsv(tx2gene_output, tx2gene_rows, PHASE7C_TX2GENE_COLUMNS)
    return tx2gene_rows, stats


def candidate_reference_status(
    phase7b_status: str,
    quantification_gene_ids: list[str],
    matched_rows: list[dict[str, str]],
) -> tuple[str, str, str]:
    if not quantification_gene_ids:
        return (
            "NO_GENE_MODEL_IDS_FROM_PHASE7B",
            "high_reference_mapping_uncertainty",
            "NOT_QUANTIFIABLE_NO_GENE_MODEL_MAP",
        )
    if not matched_rows:
        return (
            "NO_TRANSCRIPT_REFERENCE_MATCH_FOR_PHASE7B_GENE_IDS",
            "high_reference_mapping_uncertainty",
            "NOT_QUANTIFIABLE_NO_TRANSCRIPT_REFERENCE_MATCH",
        )
    if phase7b_status.startswith("NOT_READY"):
        return (
            "TRANSCRIPT_REFERENCE_PRESENT_BUT_PHASE7B_NOT_READY",
            "high_locus_or_paralog_ambiguity",
            "REFERENCE_PRESENT_BUT_NOT_READY_FOR_EXPRESSION_INTERPRETATION",
        )
    if phase7b_status.startswith("TARGETED"):
        return (
            "TRANSCRIPT_REFERENCE_PRESENT_TARGETED_REVIEW_REQUIRED",
            "high_locus_or_paralog_ambiguity",
            "REFERENCE_PRESENT_TARGETED_LOCUS_VALIDATION_REQUIRED",
        )
    matched_gene_count = len({row["gene_id"] for row in matched_rows})
    if len(matched_rows) > matched_gene_count:
        return (
            "TRANSCRIPT_REFERENCE_PRESENT_WITH_ISOFORMS",
            "moderate_isoform_summarization_risk",
            "REFERENCE_PRESENT_WITH_ISOFORM_SUMMARIZATION_CAVEAT",
        )
    return (
        "TRANSCRIPT_REFERENCE_PRESENT",
        "low_reference_mapping_risk",
        "REFERENCE_READY_FOR_FIRST_PASS_QUANTIFICATION",
    )


def build_candidate_reference_validation(
    candidate_quant_map: Path,
    tx2gene_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    candidate_rows = read_tsv(candidate_quant_map, PHASE7B_CANDIDATE_QUANT_MAP_COLUMNS)
    transcripts_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tx2gene_rows:
        transcripts_by_gene[row["gene_id"]].append(row)

    rows: list[dict[str, str]] = []
    for candidate in candidate_rows:
        gene_ids = split_delimited_values(candidate["quantification_gene_ids"])
        matched_rows: list[dict[str, str]] = []
        for gene_id in gene_ids:
            matched_rows.extend(transcripts_by_gene.get(gene_id, []))
        status, ambiguity, readiness = candidate_reference_status(
            candidate["quantification_strategy_status"],
            gene_ids,
            matched_rows,
        )
        rows.append(
            {
                "mechanism": candidate["mechanism"],
                "gene_symbol": candidate["gene_symbol"],
                "phase7b_quantification_strategy_status": candidate["quantification_strategy_status"],
                "quantification_gene_ids": join_values(gene_ids),
                "matched_gene_ids": join_values(row["gene_id"] for row in matched_rows),
                "matched_transcript_ids": join_values(row["transcript_id"] for row in matched_rows),
                "matched_transcript_count": str(len(matched_rows)),
                "total_matched_transcript_bp": str(sum(int(row["transcript_length_bp"]) for row in matched_rows)),
                "reference_mapping_status": status,
                "ambiguity_status": ambiguity,
                "phase7c_quantification_readiness": readiness,
                "required_validation": "Run read QC, quantify with recorded parameters, review multi-mapping and isoform/paralog ambiguity, and keep interpretation retina-specific.",
                "claim_language_guardrail": "Reference construction does not show expression. Do not infer activation, differential expression, pathway activity, causation, functional advantage, clinical relevance, or longevity mechanism.",
                "supporting_files": "data/interim/rnaseq/phase7c/smic_tokyo.annotation_transcripts.fna;data/interim/rnaseq/phase7c/smic_tokyo.tx2gene.tsv;results/rnaseq/phase7b_candidate_quantification_map.tsv",
            }
        )
    return rows


def build_reference_qc_row(
    config: dict,
    tx2gene_rows: list[dict[str, str]],
    stats: dict[str, int],
    candidate_row_count: int,
) -> dict[str, str]:
    phase7c = config["phase7c_expression_reference_validation"]
    lengths = [int(row["transcript_length_bp"]) for row in tx2gene_rows]
    gene_count = len({row["gene_id"] for row in tx2gene_rows})
    warning = stats["missing_scaffold_count"] > 0 or stats["zero_length_transcript_count"] > 0
    return {
        "reference_id": "SMIC_TOKYO_PHASE7C_ANNOTATION_TRANSCRIPTS",
        "annotation_gff": phase7c["annotation_gff"],
        "genome_fasta": phase7c["genome_fasta"],
        "transcript_fasta": phase7c["transcript_fasta_output"],
        "tx2gene": phase7c["tx2gene_output"],
        "transcript_count": str(len(tx2gene_rows)),
        "gene_count": str(gene_count),
        "exon_derived_transcript_count": str(stats["exon_derived_transcript_count"]),
        "cds_fallback_transcript_count": str(stats["cds_fallback_transcript_count"]),
        "missing_scaffold_count": str(stats["missing_scaffold_count"]),
        "zero_length_transcript_count": str(stats["zero_length_transcript_count"]),
        "total_transcript_bp": str(stats["total_transcript_bp"]),
        "median_transcript_length_bp": str(int(median(lengths))),
        "candidate_quant_map_rows": str(candidate_row_count),
        "validation_status": "REFERENCE_CONSTRUCTED_WITH_WARNINGS" if warning else "REFERENCE_CONSTRUCTED_AND_INPUTS_VALIDATED",
        "construction_method": phase7c.get("construction_method", "Pure-Python GFF3 exon extraction from genome FASTA."),
        "limitations": "Transcript sequences are reconstructed from the supplied annotation and assembly. This does not validate transcript models, expression, isoform usage, or locus-specific mappability.",
        "claim_language_guardrail": "Do not treat transcript-reference presence as expression evidence, pathway activity, causation, or longevity-mechanism support.",
    }


def write_report(
    qc_row: dict[str, str],
    candidate_rows: list[dict[str, str]],
    output: Path,
) -> None:
    ready = [row for row in candidate_rows if row["phase7c_quantification_readiness"] == "REFERENCE_READY_FOR_FIRST_PASS_QUANTIFICATION"]
    caveat = [row for row in candidate_rows if row["phase7c_quantification_readiness"].startswith("REFERENCE_PRESENT_WITH")]
    targeted = [row for row in candidate_rows if row["phase7c_quantification_readiness"] == "REFERENCE_PRESENT_TARGETED_LOCUS_VALIDATION_REQUIRED"]
    not_ready = [row for row in candidate_rows if row["phase7c_quantification_readiness"].startswith("NOT_QUANTIFIABLE") or row["phase7c_quantification_readiness"].startswith("REFERENCE_PRESENT_BUT_NOT_READY")]
    lines = [
        "# Phase 7c Expression Reference Construction",
        "",
        "Phase 7c constructs an annotation-derived transcript FASTA and tx2gene table, then validates candidate-panel gene IDs against that reference. It does not download reads or quantify expression.",
        "",
        "## Reference QC",
        "",
        f"- Transcript sequences: {qc_row['transcript_count']}",
        f"- Genes represented: {qc_row['gene_count']}",
        f"- Exon-derived transcripts: {qc_row['exon_derived_transcript_count']}",
        f"- CDS-fallback transcripts: {qc_row['cds_fallback_transcript_count']}",
        f"- Missing scaffold count: {qc_row['missing_scaffold_count']}",
        f"- Validation status: {qc_row['validation_status']}",
        "",
        "## Candidate Reference Readiness",
        "",
        f"- Ready for first-pass quantification: {len(ready)}",
        f"- Ready with isoform caveat: {len(caveat)}",
        f"- Targeted locus validation required: {len(targeted)}",
        f"- Not quantifiable or not interpretable yet: {len(not_ready)}",
        "",
        "Guardrail: these outputs validate reference construction and candidate ID mapping only. They do not show that any gene is detected or expressed.",
        "",
        "Supporting tables:",
        "",
        "- `data/interim/rnaseq/phase7c/smic_tokyo.annotation_transcripts.fna`",
        "- `data/interim/rnaseq/phase7c/smic_tokyo.tx2gene.tsv`",
        "- `results/rnaseq/phase7c_expression_reference_qc.tsv`",
        "- `results/rnaseq/phase7c_candidate_reference_validation.tsv`",
    ]
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase7c(
    config_path: Path,
    transcript_fasta_output: Path,
    tx2gene_output: Path,
    reference_qc_output: Path,
    candidate_validation_output: Path,
    report_output: Path,
) -> None:
    config = read_yaml(config_path)
    if "phase7c_expression_reference_validation" not in config:
        raise ValueError("config.yaml is missing phase7c_expression_reference_validation")
    phase7c = config["phase7c_expression_reference_validation"]
    annotation_gff = Path(phase7c["annotation_gff"])
    genome_fasta = Path(phase7c["genome_fasta"])
    candidate_quant_map = Path(phase7c["phase7b_candidate_quant_map"])
    for required in [annotation_gff, genome_fasta, candidate_quant_map]:
        if not required.exists():
            raise FileNotFoundError(required)

    transcript_types = set(phase7c.get("transcript_feature_types", ["mRNA", "transcript", "lncRNA"]))
    interval_feature = phase7c.get("interval_feature_preference", "exon")
    fallback_feature = phase7c.get("fallback_interval_feature", "CDS")
    LOGGER.info("Parsing transcript models from %s", annotation_gff)
    models = parse_gff_models(annotation_gff, transcript_types, interval_feature, fallback_feature)
    LOGGER.info("Constructing transcript FASTA from %d transcript models", len(models))
    tx2gene_rows, stats = write_transcript_reference(genome_fasta, transcript_fasta_output, tx2gene_output, models)
    candidate_rows = build_candidate_reference_validation(candidate_quant_map, tx2gene_rows)
    qc_row = build_reference_qc_row(config, tx2gene_rows, stats, len(candidate_rows))
    write_tsv(reference_qc_output, [qc_row], PHASE7C_REFERENCE_QC_COLUMNS)
    write_tsv(candidate_validation_output, candidate_rows, PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS)
    write_report(qc_row, candidate_rows, report_output)
    LOGGER.info("Wrote Phase 7c reference outputs with %d transcripts and %d candidate rows", len(tx2gene_rows), len(candidate_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Construct and validate the Phase 7c expression reference.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--transcript-fasta-output", type=Path, required=True)
    parser.add_argument("--tx2gene-output", type=Path, required=True)
    parser.add_argument("--reference-qc-output", type=Path, required=True)
    parser.add_argument("--candidate-validation-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_phase7c(
        args.config,
        args.transcript_fasta_output,
        args.tx2gene_output,
        args.reference_qc_output,
        args.candidate_validation_output,
        args.report_output,
    )


if __name__ == "__main__":
    main()
