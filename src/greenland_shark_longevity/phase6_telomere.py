"""Phase 6 telomere motif and telomere-gene readiness layer."""

from __future__ import annotations

import argparse
import gzip
import logging
from pathlib import Path

from .evidence import validate_guardrail_language
from .schemas import (
    COPY_NUMBER_COLUMNS,
    DUPLICATION_AUDIT_COLUMNS,
    EVIDENCE_COLUMNS,
    PHASE6_SCAFFOLD_END_ENRICHMENT_COLUMNS,
    PHASE6_TELOMERE_GENE_AUDIT_COLUMNS,
    PHASE6_TELOMERE_MOTIF_SCAN_COLUMNS,
)
from .utils import ensure_parent, read_tsv, read_yaml, write_tsv

LOGGER = logging.getLogger(__name__)
NOT_ASSESSED = "NOT_ASSESSED"


def parse_int(value: str | None, default: int = 0) -> int:
    if value in {None, "", NOT_ASSESSED}:
        return default
    return int(float(str(value)))


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def iter_fasta_terminal_context(path: Path, context_bp: int):
    header: str | None = None
    prefix = ""
    suffix = ""
    length = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, length, prefix.upper(), suffix.upper()
                header = line[1:].strip().split()[0]
                prefix = ""
                suffix = ""
                length = 0
            else:
                chunk = line.strip()
                if not chunk:
                    continue
                length += len(chunk)
                if len(prefix) < context_bp:
                    needed = context_bp - len(prefix)
                    prefix += chunk[:needed]
                suffix = (suffix + chunk)[-context_bp:]
        if header is not None:
            yield header, length, prefix.upper(), suffix.upper()


def count_overlapping(sequence: str, motif: str) -> int:
    if not sequence or not motif:
        return 0
    count = 0
    start = 0
    while True:
        idx = sequence.find(motif, start)
        if idx == -1:
            return count
        count += 1
        start = idx + 1


def motif_counts(sequence: str, motifs: list[str]) -> dict[str, int]:
    return {motif: count_overlapping(sequence, motif) for motif in motifs}


def density_per_kb(count: int, bp: int) -> float:
    if bp <= 0:
        return 0.0
    return count / (bp / 1000.0)


def enrichment_ratio(terminal_density: float, control_density: float) -> str:
    if control_density == 0.0:
        return "INF" if terminal_density > 0 else "NOT_ASSESSED"
    return f"{terminal_density / control_density:.6g}"


def build_motif_scan_rows(
    genome_fasta: Path,
    resource_id: str,
    motifs: list[str],
    terminal_window_bp: int,
    min_sequence_length_bp: int,
) -> list[dict[str, str]]:
    if not genome_fasta.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {genome_fasta}")
    if "TTAGGG" not in motifs or "CCCTAA" not in motifs:
        raise ValueError("Phase 6 MVP expects motifs to include TTAGGG and CCCTAA")

    rows: list[dict[str, str]] = []
    context_bp = terminal_window_bp * 2
    for sequence_id, length, prefix, suffix in iter_fasta_terminal_context(genome_fasta, context_bp):
        if length < min_sequence_length_bp:
            continue
        effective_window = min(terminal_window_bp, max(1, length // 2))
        left_terminal = prefix[:effective_window]
        right_terminal = suffix[-effective_window:]
        left_control = prefix[effective_window : effective_window * 2]
        suffix_genome_start = length - len(suffix)
        right_control_start = max(0, length - effective_window * 2)
        right_control_end = max(0, length - effective_window)
        right_control = suffix[
            max(0, right_control_start - suffix_genome_start) : max(0, right_control_end - suffix_genome_start)
        ]

        left_counts = motif_counts(left_terminal, motifs)
        right_counts = motif_counts(right_terminal, motifs)
        left_control_counts = motif_counts(left_control, motifs)
        right_control_counts = motif_counts(right_control, motifs)
        terminal_total = sum(left_counts.values()) + sum(right_counts.values())
        control_total = sum(left_control_counts.values()) + sum(right_control_counts.values())
        terminal_bp = len(left_terminal) + len(right_terminal)
        control_bp = len(left_control) + len(right_control)
        terminal_density = density_per_kb(terminal_total, terminal_bp)
        control_density = density_per_kb(control_total, control_bp)
        rows.append(
            {
                "resource_id": resource_id,
                "genome_fasta": str(genome_fasta),
                "sequence_id": sequence_id,
                "sequence_length_bp": str(length),
                "terminal_window_bp": str(terminal_window_bp),
                "effective_terminal_window_bp": str(effective_window),
                "left_TTAGGG_count": str(left_counts.get("TTAGGG", 0)),
                "left_CCCTAA_count": str(left_counts.get("CCCTAA", 0)),
                "right_TTAGGG_count": str(right_counts.get("TTAGGG", 0)),
                "right_CCCTAA_count": str(right_counts.get("CCCTAA", 0)),
                "left_control_TTAGGG_count": str(left_control_counts.get("TTAGGG", 0)),
                "left_control_CCCTAA_count": str(left_control_counts.get("CCCTAA", 0)),
                "right_control_TTAGGG_count": str(right_control_counts.get("TTAGGG", 0)),
                "right_control_CCCTAA_count": str(right_control_counts.get("CCCTAA", 0)),
                "terminal_total_motif_count": str(terminal_total),
                "control_total_motif_count": str(control_total),
                "terminal_motif_density_per_kb": f"{terminal_density:.6g}",
                "control_motif_density_per_kb": f"{control_density:.6g}",
                "end_enrichment_ratio": enrichment_ratio(terminal_density, control_density),
                "motif_scan_status": "CANONICAL_MOTIF_COUNTS_RECORDED_SEQUENCE_CONTEXT_ONLY",
                "interpretation_guardrail": "Do not interpret exact telomeric motif counts as telomere length, telomerase activity, rejuvenation, or longevity mechanism.",
            }
        )
    LOGGER.info("Scanned telomeric motifs in %d sequences from %s", len(rows), genome_fasta)
    return rows


def summarize_enrichment(rows: list[dict[str, str]], resource_id: str, genome_fasta: Path, motifs: list[str], terminal_window_bp: int) -> list[dict[str, str]]:
    terminal_total = sum(parse_int(row["terminal_total_motif_count"]) for row in rows)
    control_total = sum(parse_int(row["control_total_motif_count"]) for row in rows)
    terminal_bp = sum(parse_int(row["effective_terminal_window_bp"]) * 2 for row in rows)
    control_bp = terminal_bp
    terminal_density = density_per_kb(terminal_total, terminal_bp)
    control_density = density_per_kb(control_total, control_bp)
    with_terminal = sum(1 for row in rows if parse_int(row["terminal_total_motif_count"]) > 0)
    enriched = 0
    for row in rows:
        ratio = row["end_enrichment_ratio"]
        if ratio == "INF" or (ratio not in {"", NOT_ASSESSED} and float(ratio) > 1.0):
            enriched += 1
    summary = {
        "resource_id": resource_id,
        "genome_fasta": str(genome_fasta),
        "motifs": ";".join(motifs),
        "terminal_window_bp": str(terminal_window_bp),
        "sequence_count_scanned": str(len(rows)),
        "total_terminal_bp_scanned": str(terminal_bp),
        "total_control_bp_scanned": str(control_bp),
        "terminal_total_motif_count": str(terminal_total),
        "control_total_motif_count": str(control_total),
        "terminal_motif_density_per_kb": f"{terminal_density:.6g}",
        "control_motif_density_per_kb": f"{control_density:.6g}",
        "end_enrichment_ratio": enrichment_ratio(terminal_density, control_density),
        "scaffolds_with_terminal_motif": str(with_terminal),
        "scaffolds_with_terminal_motif_fraction": f"{(with_terminal / len(rows)) if rows else 0.0:.6g}",
        "scaffolds_with_terminal_enrichment_gt_1": str(enriched),
        "motif_scan_status": "CANONICAL_MOTIF_END_ENRICHMENT_RECORDED_SEQUENCE_CONTEXT_ONLY",
        "biological_interpretation": "Canonical telomeric motif counts were recorded at scaffold ends as sequence-context evidence only.",
        "claim_language_guardrail": "Do not infer telomere length, telomerase activity, rejuvenation, pathway activity, causation, or longevity mechanism from Phase 6 motif counts.",
        "required_validation": "Validate scaffold-end behavior with independent assemblies, telomere-specific methods, and manual inspection before telomere biology interpretation.",
        "supporting_files": "results/telomere/phase6_telomeric_motif_scan.tsv",
    }
    validate_guardrail_language(
        [
            {
                "gene_or_pathway": "telomeric_motif_scan",
                "biological_interpretation": summary["biological_interpretation"],
                "relevance_to_aging_longevity": NOT_ASSESSED,
                "translational_relevance": NOT_ASSESSED,
            }
        ]
    )
    return [summary]


def load_telomere_panel(candidate_panel_path: Path) -> list[dict[str, str]]:
    panel = read_yaml(candidate_panel_path)
    telomere = panel.get("panels", {}).get("telomere_shelterin")
    if not telomere:
        raise ValueError(f"{candidate_panel_path} does not define panels.telomere_shelterin")
    candidates = telomere.get("candidates", [])
    rows: list[dict[str, str]] = []
    for candidate in candidates:
        rows.append(
            {
                "mechanism": telomere["mechanism"],
                "gene_symbol": candidate["gene_symbol"],
                "synonyms": ";".join(candidate.get("synonyms", [])) if candidate.get("synonyms") else NOT_ASSESSED,
                "panel_status": "CURATED_PHASE6_TELOMERE_PANEL",
            }
        )
    return rows


def index_by_gene(path: Path, required_columns: list[str], gene_column: str) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    return {row[gene_column]: row for row in read_tsv(path, required_columns)}


def readiness_status(copy_row: dict[str, str] | None, evidence_row: dict[str, str] | None) -> str:
    if evidence_row is None and copy_row is None:
        return "PANEL_ONLY_NOT_ASSESSED_IN_CURRENT_ORTHOLOGY_OUTPUT"
    if copy_row and copy_row.get("mapping_status", "").startswith("ANNOTATION_UNCERTAINTY"):
        return "ANNOTATION_UNCERTAINTY_REQUIRES_RESCUE_OR_MANUAL_REVIEW"
    if evidence_row and evidence_row.get("evidence_tier") in {"Tier 1", "Tier 2"}:
        return "ORTHOLOGY_EVIDENCE_READY_FOR_DOMAIN_AND_CONTEXT_VALIDATION"
    if evidence_row:
        return "CURRENT_EVIDENCE_ARTIFACT_OR_INCOMPLETE"
    return "COPY_NUMBER_TABLE_ONLY_REQUIRES_EVIDENCE_INTEGRATION"


def interpretation_for_readiness(status: str) -> str:
    if status == "ORTHOLOGY_EVIDENCE_READY_FOR_DOMAIN_AND_CONTEXT_VALIDATION":
        return "Telomere-related candidate has current orthology/evidence-table support, but this does not indicate telomere length or telomerase activity."
    if status == "ANNOTATION_UNCERTAINTY_REQUIRES_RESCUE_OR_MANUAL_REVIEW":
        return "Current annotation/evidence tables do not resolve this telomere-related candidate; do not infer gene absence."
    if status == "CURRENT_EVIDENCE_ARTIFACT_OR_INCOMPLETE":
        return "Current evidence is incomplete or artifact-prone and requires validation before interpretation."
    return "Candidate is curated in the panel but has not yet been assessed in the current orthology/evidence outputs."


def build_telomere_gene_audit(
    candidate_panel_path: Path,
    copy_number_path: Path,
    duplication_audit_path: Path,
    integrated_evidence_path: Path,
) -> list[dict[str, str]]:
    panel_rows = load_telomere_panel(candidate_panel_path)
    copy_by_gene = index_by_gene(copy_number_path, COPY_NUMBER_COLUMNS, "gene_symbol")
    audit_by_gene = index_by_gene(duplication_audit_path, DUPLICATION_AUDIT_COLUMNS, "gene_symbol")
    evidence_by_gene = index_by_gene(integrated_evidence_path, EVIDENCE_COLUMNS, "gene_or_pathway")
    rows: list[dict[str, str]] = []
    for panel_row in panel_rows:
        gene = panel_row["gene_symbol"]
        copy_row = copy_by_gene.get(gene)
        audit_row = audit_by_gene.get(gene)
        evidence_row = evidence_by_gene.get(gene)
        status = readiness_status(copy_row, evidence_row)
        rows.append(
            {
                "mechanism": panel_row["mechanism"],
                "gene_symbol": gene,
                "synonyms": panel_row["synonyms"],
                "panel_status": panel_row["panel_status"],
                "copy_number_mapping_status": copy_row.get("mapping_status", NOT_ASSESSED) if copy_row else NOT_ASSESSED,
                "copy_count": copy_row.get("copy_count", NOT_ASSESSED) if copy_row else NOT_ASSESSED,
                "orthogroup_id": copy_row.get("orthogroup_id", NOT_ASSESSED) if copy_row else NOT_ASSESSED,
                "resources_supporting": evidence_row.get("resources_supporting", copy_row.get("resources_supporting", NOT_ASSESSED) if copy_row else NOT_ASSESSED) if evidence_row else (copy_row.get("resources_supporting", NOT_ASSESSED) if copy_row else NOT_ASSESSED),
                "duplication_artifact_risk": audit_row.get("artifact_risk", NOT_ASSESSED) if audit_row else NOT_ASSESSED,
                "integrated_evidence_tier": evidence_row.get("evidence_tier", NOT_ASSESSED) if evidence_row else NOT_ASSESSED,
                "integrated_artifact_risk": evidence_row.get("artifact_risk", NOT_ASSESSED) if evidence_row else NOT_ASSESSED,
                "readiness_status": status,
                "conservative_interpretation": interpretation_for_readiness(status),
                "required_validation": "Confirm orthology, domain architecture, locus context, annotation consistency, and cross-resource support before telomere-related biological interpretation.",
                "claim_language_guardrail": "Do not infer telomere length, telomerase activity, rejuvenation, gene absence, pathway activity, causation, or longevity mechanism from Phase 6 gene readiness.",
                "supporting_files": "config/candidate_panels.yaml;results/orthology/candidate_copy_number.tsv;results/evidence/integrated_evidence.tsv",
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


def write_report(
    enrichment_rows: list[dict[str, str]],
    gene_rows: list[dict[str, str]],
    output: Path,
) -> None:
    summary = enrichment_rows[0]
    lines = [
        "# Phase 6 Telomere Readiness Report",
        "",
        "Phase 6 records exact canonical telomeric motifs at scaffold ends and audits telomere/shelterin candidate readiness. It does not infer telomere length, telomerase activity, rejuvenation, or longevity mechanism.",
        "",
        "## Motif Scan",
        "",
        f"- Sequences scanned: {summary['sequence_count_scanned']}",
        f"- Terminal window: {summary['terminal_window_bp']} bp",
        f"- Motifs: `{summary['motifs']}`",
        f"- Terminal motif density per kb: {summary['terminal_motif_density_per_kb']}",
        f"- Control motif density per kb: {summary['control_motif_density_per_kb']}",
        f"- End enrichment ratio: {summary['end_enrichment_ratio']}",
        f"- Scaffolds with at least one terminal motif: {summary['scaffolds_with_terminal_motif']} ({summary['scaffolds_with_terminal_motif_fraction']})",
        "",
        "## Telomere-Gene Readiness",
        "",
        "| Gene | Current tier | Readiness | Artifact risk |",
        "|---|---|---|---|",
    ]
    for row in gene_rows:
        lines.append(
            f"| `{row['gene_symbol']}` | {row['integrated_evidence_tier']} | {row['readiness_status']} | {row['integrated_artifact_risk']} |"
        )
    lines.extend(
        [
            "",
            "Guardrail: motif enrichment and candidate-gene readiness are sequence/context observations only.",
            "",
            "Supporting tables:",
            "",
            "- `results/telomere/phase6_telomeric_motif_scan.tsv`",
            "- `results/telomere/phase6_scaffold_end_enrichment.tsv`",
            "- `results/telomere/phase6_telomere_gene_audit.tsv`",
        ]
    )
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_phase6_outputs(
    config_path: Path,
    motif_scan_output: Path,
    enrichment_output: Path,
    gene_audit_output: Path,
    report_output: Path,
) -> None:
    config = read_yaml(config_path)
    phase6 = config.get("phase6_telomere", {})
    genome_fasta = Path(phase6.get("genome_fasta", ""))
    resource_id = str(phase6.get("resource_id", "SMIC_TOKYO_GENOME_2025"))
    motifs = [str(motif).upper() for motif in phase6.get("motifs", ["TTAGGG", "CCCTAA"])]
    terminal_window_bp = parse_int(str(phase6.get("terminal_window_bp", 10000)), default=10000)
    min_sequence_length_bp = parse_int(str(phase6.get("min_sequence_length_bp", 1000)), default=1000)

    motif_rows = build_motif_scan_rows(genome_fasta, resource_id, motifs, terminal_window_bp, min_sequence_length_bp)
    enrichment_rows = summarize_enrichment(motif_rows, resource_id, genome_fasta, motifs, terminal_window_bp)
    gene_rows = build_telomere_gene_audit(
        Path(phase6.get("candidate_panel", "config/candidate_panels.yaml")),
        Path(phase6.get("candidate_copy_number", "results/orthology/candidate_copy_number.tsv")),
        Path(phase6.get("candidate_duplication_audit", "results/orthology/candidate_duplication_audit.tsv")),
        Path(phase6.get("integrated_evidence", "results/evidence/integrated_evidence.tsv")),
    )

    write_tsv(motif_scan_output, motif_rows, PHASE6_TELOMERE_MOTIF_SCAN_COLUMNS)
    write_tsv(enrichment_output, enrichment_rows, PHASE6_SCAFFOLD_END_ENRICHMENT_COLUMNS)
    write_tsv(gene_audit_output, gene_rows, PHASE6_TELOMERE_GENE_AUDIT_COLUMNS)
    write_report(enrichment_rows, gene_rows, report_output)
    LOGGER.info("Wrote Phase 6 telomere outputs for %d sequences and %d genes", len(motif_rows), len(gene_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 6 telomere motif and gene-readiness layer.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--motif-scan-output", type=Path, required=True)
    parser.add_argument("--enrichment-output", type=Path, required=True)
    parser.add_argument("--gene-audit-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    write_phase6_outputs(
        args.config,
        args.motif_scan_output,
        args.enrichment_output,
        args.gene_audit_output,
        args.report_output,
    )


if __name__ == "__main__":
    main()
