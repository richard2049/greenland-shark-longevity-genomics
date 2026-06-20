"""Phase 4b genome-level validation for targeted annotation rescue candidates."""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .schemas import (
    PHASE4_RESCUE_QUERY_COLUMNS,
    PHASE4_RESCUE_SUMMARY_COLUMNS,
    PHASE4B_ALIGNMENT_HIT_COLUMNS,
    PHASE4B_PREFLIGHT_COLUMNS,
    PHASE4B_QUERY_COLUMNS,
    PHASE4B_SUMMARY_COLUMNS,
    PHASE4B_TARGET_REGION_COLUMNS,
)
from .utils import configure_logging, read_tsv, read_yaml, write_tsv


PROTEIN_RESCUE_CANDIDATE = "PROTEIN_LEVEL_RESCUE_CANDIDATE"
NOT_RUN = "NOT_RUN_REQUIRES_INPUTS"


@dataclass(frozen=True)
class ProteinRecord:
    protein_id: str
    description: str
    sequence: str
    source_fasta: str

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class QueryRecord:
    gene_symbol: str
    mechanism: str
    query_id: str
    query_source_type: str
    original_protein_id: str
    original_description: str
    source_species_id: str
    source_resource_id: str
    source_fasta: str
    sequence: str
    selection_reason: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    def row(self) -> dict[str, str]:
        return {
            "gene_symbol": self.gene_symbol,
            "mechanism": self.mechanism,
            "query_id": self.query_id,
            "query_source_type": self.query_source_type,
            "original_protein_id": self.original_protein_id,
            "original_description": self.original_description,
            "source_species_id": self.source_species_id,
            "source_resource_id": self.source_resource_id,
            "source_fasta": self.source_fasta,
            "sequence_length_aa": str(self.length),
            "selection_reason": self.selection_reason,
            "query_status": "READY_FOR_GENOME_ALIGNMENT" if self.sequence else "SEQUENCE_UNAVAILABLE",
            "notes": "Phase 4b query; alignment result is required before genome-level interpretation.",
        }


@dataclass(frozen=True)
class AlignmentHit:
    gene_symbol: str
    query_id: str
    query_source_type: str
    original_protein_id: str
    source_species_id: str
    source_resource_id: str
    target_seqid: str
    target_start: int
    target_end: int
    strand: str
    feature_id: str
    feature_type: str
    query_start: int
    query_end: int
    query_length_aa: int
    query_coverage: float
    cds_feature_count: int
    exon_feature_count: int
    frameshift_or_stop_flag: bool
    identity: str
    positive: str
    score: str
    rank: str
    raw_attributes: str

    @property
    def locus(self) -> tuple[str, int, int]:
        return (self.target_seqid, min(self.target_start, self.target_end), max(self.target_start, self.target_end))


def parse_fasta(path: Path) -> dict[str, ProteinRecord]:
    records: dict[str, ProteinRecord] = {}
    current_header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    record = record_from_header(current_header, chunks, path)
                    records[record.protein_id] = record
                current_header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if current_header is not None:
        record = record_from_header(current_header, chunks, path)
        records[record.protein_id] = record
    return records


def record_from_header(header: str, chunks: list[str], path: Path) -> ProteinRecord:
    protein_id, _, description = header.partition(" ")
    return ProteinRecord(protein_id, description, sanitize_sequence("".join(chunks)), str(path))


def sanitize_sequence(sequence: str) -> str:
    return "".join(residue if "A" <= residue <= "Z" else "X" for residue in sequence.upper().replace("*", ""))


def split_assessed(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip() and part.strip() != "NOT_ASSESSED"]


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.:-]+", "_", value.strip())
    return token.strip("_") or "UNKNOWN"


def make_query_id(gene_symbol: str, source_type: str, protein_id: str) -> str:
    return f"{safe_token(gene_symbol)}__{safe_token(source_type)}__{safe_token(protein_id)}"


def load_staged_protein(record_id: str, source_fasta: Path) -> ProteinRecord | None:
    records = parse_fasta(source_fasta)
    if record_id in records:
        return records[record_id]
    normalized = record_id.replace("|", "_").replace(":", "_")
    for protein_id, record in records.items():
        if protein_id.replace("|", "_").replace(":", "_") == normalized:
            return record
    return None


def build_query_records(
    summary_path: Path,
    query_inventory_path: Path,
    rescue_candidate_fasta: Path,
    include_comparator_for_rescued: bool = False,
) -> list[QueryRecord]:
    summary_rows = read_tsv(summary_path, PHASE4_RESCUE_SUMMARY_COLUMNS)
    selected_query_rows = [
        row
        for row in read_tsv(query_inventory_path, PHASE4_RESCUE_QUERY_COLUMNS)
        if row["selection_status"] == "SELECTED_FOR_RESCUE_QUERY"
    ]
    rescue_sequences = parse_fasta(rescue_candidate_fasta) if rescue_candidate_fasta.exists() else {}
    queries: list[QueryRecord] = []
    seen: set[str] = set()

    for row in summary_rows:
        gene_symbol = row["gene_symbol"]
        if row["protein_rescue_status"] == PROTEIN_RESCUE_CANDIDATE:
            for protein_id in split_assessed(row["candidate_protein_ids"]):
                record = rescue_sequences.get(protein_id)
                if record is None:
                    continue
                query_id = make_query_id(gene_symbol, "focal_rescue_candidate", protein_id)
                if query_id in seen:
                    continue
                seen.add(query_id)
                queries.append(
                    QueryRecord(
                        gene_symbol=gene_symbol,
                        mechanism=row["mechanism"],
                        query_id=query_id,
                        query_source_type="focal_rescue_candidate",
                        original_protein_id=protein_id,
                        original_description=record.description,
                        source_species_id="smic",
                        source_resource_id="SMIC_TOKYO_GENOME_2025",
                        source_fasta=record.source_fasta,
                        sequence=record.sequence,
                        selection_reason="reciprocal protein-level rescue candidate from Phase 4",
                    )
                )

        if row["protein_rescue_status"] != PROTEIN_RESCUE_CANDIDATE or include_comparator_for_rescued:
            for query_row in selected_query_rows:
                if query_row["gene_symbol"] != gene_symbol:
                    continue
                source_fasta = Path(query_row["source_fasta"])
                if not source_fasta.exists():
                    continue
                record = load_staged_protein(query_row["query_protein_id"], source_fasta)
                if record is None:
                    continue
                query_id = make_query_id(gene_symbol, "comparator_query", query_row["query_protein_id"])
                if query_id in seen:
                    continue
                seen.add(query_id)
                queries.append(
                    QueryRecord(
                        gene_symbol=gene_symbol,
                        mechanism=row["mechanism"],
                        query_id=query_id,
                        query_source_type="comparator_query_for_unresolved_gene",
                        original_protein_id=query_row["query_protein_id"],
                        original_description=query_row["query_description"],
                        source_species_id=query_row["query_species_id"],
                        source_resource_id=query_row["query_resource_id"],
                        source_fasta=query_row["source_fasta"],
                        sequence=record.sequence,
                        selection_reason="selected comparator query for unresolved Phase 4 target",
                    )
                )
    return queries


def write_query_fasta(path: Path, queries: list[QueryRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for query in sorted(queries, key=lambda item: item.query_id):
            handle.write(f">{query.query_id} {query.original_protein_id} {query.original_description}\n")
            for idx in range(0, len(query.sequence), 80):
                handle.write(query.sequence[idx : idx + 80] + "\n")


def parse_candidate_scaffolds(summary_rows: list[dict[str, str]]) -> dict[str, set[str]]:
    scaffolds_by_gene: dict[str, set[str]] = {}
    for row in summary_rows:
        gene = row["gene_symbol"]
        for part in row.get("candidate_coordinate_summary", "").split(";"):
            fields = part.split(":")
            if len(fields) >= 3 and fields[1] != "NO_COORDINATES":
                scaffolds_by_gene.setdefault(gene, set()).add(fields[1])
    return scaffolds_by_gene


def header_matches_seqid(header: str, seqid: str) -> bool:
    first_token = header.split()[0]
    if first_token == seqid:
        return True
    return re.search(rf"(^|[\s,]){re.escape(seqid)}($|[\s,])", header) is not None


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def extract_target_scaffolds(
    genome_fasta: Path,
    target_fasta: Path,
    scaffolds_by_gene: dict[str, set[str]],
) -> tuple[set[str], list[dict[str, str]]]:
    requested = sorted({seqid for seqids in scaffolds_by_gene.values() for seqid in seqids})
    requested_set = set(requested)
    found: dict[str, str] = {}
    target_fasta.parent.mkdir(parents=True, exist_ok=True)
    active_header: str | None = None
    active_seqid: str | None = None
    with open_text(genome_fasta) as source, target_fasta.open("w", encoding="utf-8") as output:
        for line in source:
            if line.startswith(">"):
                header = line[1:].strip()
                active_seqid = next((seqid for seqid in requested if header_matches_seqid(header, seqid)), None)
                active_header = header if active_seqid else None
                if active_seqid:
                    found[active_seqid] = header
                    output.write(line)
            elif active_header:
                output.write(line)

    rows: list[dict[str, str]] = []
    for gene, seqids in sorted(scaffolds_by_gene.items()):
        for seqid in sorted(seqids):
            rows.append(
                {
                    "gene_symbol": gene,
                    "source_coordinate": seqid,
                    "requested_seqid": seqid,
                    "target_fasta_header": found.get(seqid, "NOT_FOUND"),
                    "target_fasta_path": str(target_fasta),
                    "extraction_status": "EXTRACTED" if seqid in found else "NOT_FOUND",
                    "notes": "Target scaffold extracted for limited-scope genome validation.",
                }
            )
    missing = requested_set - set(found)
    for seqid in sorted(missing):
        if not any(row["requested_seqid"] == seqid for row in rows):
            rows.append(
                {
                    "gene_symbol": "NOT_ASSESSED",
                    "source_coordinate": seqid,
                    "requested_seqid": seqid,
                    "target_fasta_header": "NOT_FOUND",
                    "target_fasta_path": str(target_fasta),
                    "extraction_status": "NOT_FOUND",
                    "notes": "Requested target scaffold was not found in the genome FASTA headers.",
                }
            )
    return set(found), rows


def resolve_miniprot(executable: str) -> str | None:
    path = Path(executable)
    if path.exists():
        return str(path)
    return shutil.which(executable)


def build_preflight_rows(
    genome_fasta: Path,
    target_fasta: Path,
    genome_scope: str,
    miniprot_executable: str,
    queries: list[QueryRecord],
    raw_gff: Path,
) -> tuple[list[dict[str, str]], bool, str | None]:
    miniprot_path = resolve_miniprot(miniprot_executable)
    ready = target_fasta.exists() and bool(miniprot_path) and bool(queries)
    rows = [
        {
            "check_id": "genome_fasta",
            "resource": "target_genome",
            "path_or_command": str(genome_fasta),
            "status": "LOCAL_AVAILABLE" if genome_fasta.exists() else "MISSING_LOCAL",
            "details": str(genome_fasta.stat().st_size) if genome_fasta.exists() else "NOT_ASSESSED",
            "required_for_alignment": "true",
            "notes": "Genome FASTA is required for protein-to-genome rescue.",
        },
        {
            "check_id": "alignment_target_fasta",
            "resource": "target_genome_scope",
            "path_or_command": str(target_fasta),
            "status": "LOCAL_AVAILABLE" if target_fasta.exists() else "MISSING_LOCAL",
            "details": genome_scope,
            "required_for_alignment": "true",
            "notes": "This is the actual FASTA passed to miniprot. It may be a scaffold subset for laptop-scale validation.",
        },
        {
            "check_id": "miniprot_executable",
            "resource": "miniprot",
            "path_or_command": miniprot_executable,
            "status": "LOCAL_AVAILABLE" if miniprot_path else "MISSING_LOCAL",
            "details": miniprot_path or "NOT_ASSESSED",
            "required_for_alignment": "true",
            "notes": "miniprot is required for genome-level protein-to-genome alignment.",
        },
        {
            "check_id": "phase4b_queries",
            "resource": "query_proteins",
            "path_or_command": "data/interim/annotation_rescue/phase4b_miniprot_queries.faa",
            "status": "LOCAL_AVAILABLE" if queries else "MISSING_LOCAL",
            "details": str(len(queries)),
            "required_for_alignment": "true",
            "notes": "Queries include focal protein-level rescue candidates and comparator queries for unresolved genes.",
        },
        {
            "check_id": "raw_gff_output",
            "resource": "workflow_output",
            "path_or_command": str(raw_gff),
            "status": "READY_TO_WRITE",
            "details": "NOT_ASSESSED",
            "required_for_alignment": "false",
            "notes": "Raw miniprot GFF output path.",
        },
    ]
    return rows, ready, miniprot_path


def parse_attributes(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in value.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, attr_value = part.split("=", 1)
        elif " " in part:
            key, attr_value = part.split(" ", 1)
        else:
            key, attr_value = part, ""
        attrs[key.strip()] = attr_value.strip().strip('"')
    return attrs


def parse_target(value: str) -> tuple[str, int, int]:
    parts = value.split()
    if len(parts) < 3:
        return ("NOT_ASSESSED", 0, 0)
    return (parts[0], parse_int(parts[1]), parse_int(parts[2]))


def parse_int(value: str) -> int:
    try:
        return int(float(value))
    except ValueError:
        return 0


def disruption_flag(attrs: dict[str, str], raw_line: str, feature_type: str = "") -> bool:
    text = " ".join([raw_line, " ".join(f"{key}={value}" for key, value in attrs.items())]).lower()
    feature = feature_type.lower()
    if feature == "frameshift" or "frameshift" in text:
        return True
    return any(token in text for token in ["premature_stop", "internal_stop"])


def parse_paf_tags(line: str) -> dict[str, str]:
    fields = line.rstrip("\n").split("\t")
    tags: dict[str, str] = {}
    for field in fields[13:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def paf_disruption_flag(line: str) -> bool:
    tags = parse_paf_tags(line)
    frameshifts = parse_int(tags.get("fs", "0"))
    stops = parse_int(tags.get("st", "0"))
    return frameshifts > 0 or stops > 0


def parse_miniprot_gff(path: Path, query_lookup: dict[str, QueryRecord]) -> list[AlignmentHit]:
    parents: dict[str, dict[str, object]] = {}
    child_counts: dict[str, dict[str, int | bool]] = {}
    pending_paf_disruption = False
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##PAF\t"):
                pending_paf_disruption = paf_disruption_flag(line)
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            attrs = parse_attributes(attributes)
            feature_id = attrs.get("ID", f"{seqid}:{start}-{end}:{feature_type}")
            parent_id = attrs.get("Parent", "")
            if parent_id:
                counts = child_counts.setdefault(parent_id, {"CDS": 0, "exon": 0, "disruption": False})
                if feature_type == "CDS":
                    counts["CDS"] = int(counts["CDS"]) + 1
                if feature_type == "exon":
                    counts["exon"] = int(counts["exon"]) + 1
                counts["disruption"] = bool(counts["disruption"]) or disruption_flag(attrs, line, feature_type)

            if feature_type not in {"mRNA", "match", "protein_match"}:
                continue
            target_name, query_start, query_end = parse_target(attrs.get("Target", ""))
            query = query_lookup.get(target_name)
            if query is None:
                continue
            parents[feature_id] = {
                "gene_symbol": query.gene_symbol,
                "query_id": query.query_id,
                "query_source_type": query.query_source_type,
                "original_protein_id": query.original_protein_id,
                "source_species_id": query.source_species_id,
                "source_resource_id": query.source_resource_id,
                "target_seqid": seqid,
                "target_start": parse_int(start),
                "target_end": parse_int(end),
                "strand": strand,
                "feature_id": feature_id,
                "feature_type": feature_type,
                "query_start": query_start,
                "query_end": query_end,
                "query_length_aa": query.length,
                "identity": attrs.get("Identity", attrs.get("identity", "NOT_ASSESSED")),
                "positive": attrs.get("Positive", attrs.get("positive", "NOT_ASSESSED")),
                "score": score if score != "." else attrs.get("AS", "NOT_ASSESSED"),
                "rank": attrs.get("Rank", attrs.get("rank", "NOT_ASSESSED")),
                "raw_attributes": attributes,
                "parent_disruption": pending_paf_disruption or disruption_flag(attrs, line, feature_type),
            }
            pending_paf_disruption = False
    hits: list[AlignmentHit] = []
    for feature_id, values in parents.items():
        query_length = int(values["query_length_aa"])
        query_span = max(0, int(values["query_end"]) - int(values["query_start"]) + 1)
        coverage = query_span / query_length if query_length else 0.0
        counts = child_counts.get(feature_id, {"CDS": 0, "exon": 0, "disruption": False})
        hits.append(
            AlignmentHit(
                gene_symbol=str(values["gene_symbol"]),
                query_id=str(values["query_id"]),
                query_source_type=str(values["query_source_type"]),
                original_protein_id=str(values["original_protein_id"]),
                source_species_id=str(values["source_species_id"]),
                source_resource_id=str(values["source_resource_id"]),
                target_seqid=str(values["target_seqid"]),
                target_start=int(values["target_start"]),
                target_end=int(values["target_end"]),
                strand=str(values["strand"]),
                feature_id=str(values["feature_id"]),
                feature_type=str(values["feature_type"]),
                query_start=int(values["query_start"]),
                query_end=int(values["query_end"]),
                query_length_aa=query_length,
                query_coverage=coverage,
                cds_feature_count=int(counts["CDS"]),
                exon_feature_count=int(counts["exon"]),
                frameshift_or_stop_flag=bool(counts["disruption"]) or bool(values["parent_disruption"]),
                identity=str(values["identity"]),
                positive=str(values["positive"]),
                score=str(values["score"]),
                rank=str(values["rank"]),
                raw_attributes=str(values["raw_attributes"]),
            )
        )
    return hits


def hit_status(hit: AlignmentHit, full_coverage: float, partial_coverage: float) -> str:
    if hit.frameshift_or_stop_flag:
        return "POSSIBLE_DISRUPTION"
    if hit.query_coverage >= full_coverage:
        return "HIGH_COVERAGE_NO_DISRUPTION"
    if hit.query_coverage >= partial_coverage:
        return "PARTIAL_ALIGNMENT"
    return "LOW_COVERAGE_ALIGNMENT"


def alignment_rows(hits: list[AlignmentHit], full_coverage: float, partial_coverage: float) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for hit in hits:
        rows.append(
            {
                "gene_symbol": hit.gene_symbol,
                "query_id": hit.query_id,
                "query_source_type": hit.query_source_type,
                "original_protein_id": hit.original_protein_id,
                "source_species_id": hit.source_species_id,
                "source_resource_id": hit.source_resource_id,
                "target_seqid": hit.target_seqid,
                "target_start": str(hit.target_start),
                "target_end": str(hit.target_end),
                "strand": hit.strand,
                "feature_id": hit.feature_id,
                "feature_type": hit.feature_type,
                "query_start": str(hit.query_start),
                "query_end": str(hit.query_end),
                "query_length_aa": str(hit.query_length_aa),
                "query_coverage": f"{hit.query_coverage:.6g}",
                "cds_feature_count": str(hit.cds_feature_count),
                "exon_feature_count": str(hit.exon_feature_count),
                "frameshift_or_stop_flag": str(hit.frameshift_or_stop_flag),
                "identity": hit.identity,
                "positive": hit.positive,
                "score": hit.score,
                "rank": hit.rank,
                "alignment_status": hit_status(hit, full_coverage, partial_coverage),
                "raw_attributes": hit.raw_attributes,
                "notes": "miniprot genome alignment hit; manual locus review remains required.",
            }
        )
    return rows


def overlap_fraction(a: tuple[str, int, int], b: tuple[str, int, int]) -> float:
    if a[0] != b[0]:
        return 0.0
    left = max(a[1], b[1])
    right = min(a[2], b[2])
    if right < left:
        return 0.0
    overlap = right - left + 1
    return overlap / max(1, min(a[2] - a[1] + 1, b[2] - b[1] + 1))


def separable_loci(hits: list[AlignmentHit]) -> list[tuple[str, int, int]]:
    loci: list[tuple[str, int, int]] = []
    for hit in sorted(hits, key=lambda item: (item.target_seqid, item.target_start, item.target_end)):
        locus = hit.locus
        if all(overlap_fraction(locus, existing) <= 0.5 for existing in loci):
            loci.append(locus)
    return loci


def summary_rows(
    summary_path: Path,
    queries: list[QueryRecord],
    hits: list[AlignmentHit],
    ready: bool,
    full_coverage: float,
    partial_coverage: float,
    eligible_genes: set[str] | None = None,
    genome_scope: str = "full_genome",
) -> list[dict[str, str]]:
    phase4_rows = read_tsv(summary_path, PHASE4_RESCUE_SUMMARY_COLUMNS)
    rows: list[dict[str, str]] = []
    queries_by_gene: dict[str, list[QueryRecord]] = {}
    hits_by_gene: dict[str, list[AlignmentHit]] = {}
    for query in queries:
        queries_by_gene.setdefault(query.gene_symbol, []).append(query)
    for hit in hits:
        hits_by_gene.setdefault(hit.gene_symbol, []).append(hit)

    for row in phase4_rows:
        gene = row["gene_symbol"]
        gene_hits = hits_by_gene.get(gene, [])
        high = [
            hit
            for hit in gene_hits
            if hit.query_coverage >= full_coverage and not hit.frameshift_or_stop_flag
        ]
        partial = [hit for hit in gene_hits if hit.query_coverage >= partial_coverage]
        disrupted = [hit for hit in gene_hits if hit.frameshift_or_stop_flag]
        loci = separable_loci(high)
        best_coverage = max([hit.query_coverage for hit in gene_hits], default=0.0)
        if ready and eligible_genes is not None and gene not in eligible_genes:
            validation_status = "LIMITED_SCOPE_NO_CANDIDATE_REGION"
            classification = "not_assessable"
            artifact_risk = "not_assessable"
            notes = "Genome validation used a limited target FASTA and this gene had no protein-level rescued candidate region."
        elif not ready:
            validation_status = NOT_RUN
            classification = "not_assessable"
            artifact_risk = "not_assessable"
            notes = "Genome-level rescue was not run because required local inputs are missing."
        elif len(loci) >= 2:
            validation_status = "MULTIPLE_SEPARABLE_HIGH_COVERAGE_ALIGNMENTS"
            classification = "duplicated_candidate"
            artifact_risk = "moderate"
            notes = "Multiple separable high-coverage protein-to-genome alignments were detected; manual locus review is required."
        elif len(loci) == 1:
            validation_status = "HIGH_COVERAGE_GENOME_ALIGNMENT"
            classification = "intact_candidate"
            artifact_risk = "moderate"
            notes = "A high-coverage protein-to-genome alignment was detected without parsed disruption flags; manual exon/locus review is required."
        elif disrupted:
            validation_status = "POSSIBLE_DISRUPTED_ALIGNMENT"
            classification = "possible_pseudogene"
            artifact_risk = "high"
            notes = "Alignment output contains frameshift/stop-like flags; this is a candidate disruption signal requiring manual review."
        elif partial:
            validation_status = "PARTIAL_GENOME_ALIGNMENT"
            classification = "possible_fragment"
            artifact_risk = "high"
            notes = "Only partial protein-to-genome alignment support was detected."
        else:
            validation_status = "NO_QUALIFYING_GENOME_ALIGNMENT"
            classification = "annotation_uncertainty"
            artifact_risk = "high"
            notes = "No qualifying genome-level alignment was detected under the configured filters; do not infer absence."

        rows.append(
            {
                "gene_symbol": gene,
                "mechanism": row["mechanism"],
                "protein_rescue_status": row["protein_rescue_status"],
                "query_count": str(len(queries_by_gene.get(gene, []))),
                "alignment_hit_count": str(len(gene_hits)),
                "high_coverage_hit_count": str(len(high)),
                "separable_high_coverage_loci": str(len(loci)),
                "best_query_coverage": f"{best_coverage:.6g}" if ready else "NOT_ASSESSED",
                "candidate_loci": ";".join(f"{seqid}:{start}-{end}" for seqid, start, end in loci) if loci else "NOT_ASSESSED",
                "genome_validation_status": validation_status,
                "candidate_classification": classification,
                "artifact_risk": artifact_risk,
                "claim_language_guardrail": "Do not report absent, inactivated, activated, adaptive, causal, or validated duplication language from this table alone.",
                "required_validation": "Inspect raw miniprot GFF, exon structure, coordinates, overlap between hits, domain support, genome scope, and cross-resource evidence before biological interpretation.",
                "notes": f"{notes} Genome scope: {genome_scope}.",
            }
        )
    return rows


def miniprot_path_arg(path: Path, miniprot: str) -> str:
    # Docker wrappers mount the repository as /work, so arguments must use
    # Linux-style separators even when the workflow is launched from Windows.
    if miniprot.lower().endswith((".cmd", ".bat")):
        return path.as_posix()
    return path.as_posix()


def run_miniprot(miniprot: str, genome_fasta: Path, query_fasta: Path, raw_gff: Path, stderr_log: Path, threads: int, extra_args: list[str]) -> None:
    raw_gff.parent.mkdir(parents=True, exist_ok=True)
    stderr_log.parent.mkdir(parents=True, exist_ok=True)
    command = [
        miniprot,
        "-t",
        str(threads),
        "--gff",
        *extra_args,
        miniprot_path_arg(genome_fasta, miniprot),
        miniprot_path_arg(query_fasta, miniprot),
    ]
    logging.info("Running miniprot: %s", " ".join(command))
    with raw_gff.open("w", encoding="utf-8") as stdout, stderr_log.open("w", encoding="utf-8") as stderr:
        try:
            subprocess.run(command, check=True, stdout=stdout, stderr=stderr, text=True)
        except subprocess.CalledProcessError as error:
            stderr.flush()
            try:
                details = stderr_log.read_text(encoding="utf-8", errors="replace") if stderr_log.exists() else ""
            except OSError:
                details = "stderr log could not be read while the failing subprocess handle was open"
            if error.returncode == 137:
                hint = (
                    "miniprot exited with code 137, which usually means Docker killed the container "
                    "because of a memory limit. Increase Docker Desktop memory, reduce threads, or "
                    "increase miniprot -M sampling in config/phase4b_genome_validation.extra_miniprot_args."
                )
            else:
                hint = "miniprot failed; inspect the command, genome FASTA, query FASTA, and stderr."
            raise RuntimeError(f"{hint}\nCommand: {' '.join(command)}\nStderr:\n{details}") from error


def write_not_run_gff(path: Path, reason: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        handle.write(f"# miniprot_not_run={reason}\n")


def phase4b_config(config_path: Path) -> dict:
    config = read_yaml(config_path)
    params = config.get("phase4b_genome_validation", {})
    if not isinstance(params, dict):
        raise ValueError("config/phase4b_genome_validation must be a mapping")
    return config


def run_phase4b(
    config_path: Path,
    phase4_summary: Path,
    query_inventory: Path,
    rescue_candidate_fasta: Path,
    query_fasta_output: Path,
    target_fasta_output: Path,
    target_regions_output: Path,
    preflight_output: Path,
    raw_gff_output: Path,
    stderr_log_output: Path,
    alignment_hits_output: Path,
    summary_output: Path,
) -> None:
    config = phase4b_config(config_path)
    params = config.get("phase4b_genome_validation", {})
    genome_fasta = Path(str(params.get("genome_fasta", config.get("reference_inputs", {}).get("assembly_fasta", ""))))
    genome_scope = str(params.get("genome_scope", "candidate_scaffolds"))
    miniprot_executable = str(params.get("miniprot_executable", "miniprot"))
    include_comparator_for_rescued = bool(params.get("include_comparator_queries_for_genes_with_protein_rescue", False))
    full_coverage = float(params.get("min_query_coverage_for_intact_candidate", 0.70))
    partial_coverage = float(params.get("min_query_coverage_for_partial_candidate", 0.20))
    threads = int(params.get("threads", 1))
    fail_if_missing = bool(params.get("fail_if_missing_inputs", False))
    extra_args = [str(item) for item in params.get("extra_miniprot_args", [])]

    queries = build_query_records(
        phase4_summary,
        query_inventory,
        rescue_candidate_fasta,
        include_comparator_for_rescued=include_comparator_for_rescued,
    )
    write_query_fasta(query_fasta_output, queries)

    phase4_rows = read_tsv(phase4_summary, PHASE4_RESCUE_SUMMARY_COLUMNS)
    scaffolds_by_gene = parse_candidate_scaffolds(phase4_rows)
    eligible_genes: set[str] | None = None
    if genome_scope == "candidate_scaffolds":
        if genome_fasta.exists():
            found_scaffolds, target_region_rows = extract_target_scaffolds(genome_fasta, target_fasta_output, scaffolds_by_gene)
        else:
            found_scaffolds = set()
            target_region_rows = [
                {
                    "gene_symbol": gene,
                    "source_coordinate": seqid,
                    "requested_seqid": seqid,
                    "target_fasta_header": "NOT_ASSESSED",
                    "target_fasta_path": str(target_fasta_output),
                    "extraction_status": "GENOME_FASTA_MISSING",
                    "notes": "Target scaffold extraction requires the configured genome FASTA.",
                }
                for gene, seqids in sorted(scaffolds_by_gene.items())
                for seqid in sorted(seqids)
            ]
        eligible_genes = {
            gene
            for gene, seqids in scaffolds_by_gene.items()
            if any(seqid in found_scaffolds for seqid in seqids)
        }
        target_fasta = target_fasta_output
    elif genome_scope == "full_genome":
        target_region_rows = []
        target_fasta = genome_fasta
    else:
        raise ValueError("phase4b_genome_validation.genome_scope must be candidate_scaffolds or full_genome")
    write_tsv(target_regions_output, target_region_rows, PHASE4B_TARGET_REGION_COLUMNS)

    preflight_rows, ready, miniprot_path = build_preflight_rows(
        genome_fasta,
        target_fasta,
        genome_scope,
        miniprot_executable,
        queries,
        raw_gff_output,
    )
    write_tsv(preflight_output, preflight_rows, PHASE4B_PREFLIGHT_COLUMNS)
    write_tsv(query_fasta_output.with_suffix(".tsv"), [query.row() for query in queries], PHASE4B_QUERY_COLUMNS)

    hits: list[AlignmentHit] = []
    if ready and miniprot_path:
        run_miniprot(miniprot_path, target_fasta, query_fasta_output, raw_gff_output, stderr_log_output, threads, extra_args)
        query_lookup = {query.query_id: query for query in queries}
        hits = parse_miniprot_gff(raw_gff_output, query_lookup)
    else:
        reason = "missing genome FASTA, miniprot executable, or query proteins"
        if fail_if_missing:
            raise RuntimeError(f"Phase 4b genome validation cannot run: {reason}")
        write_not_run_gff(raw_gff_output, reason)
        stderr_log_output.parent.mkdir(parents=True, exist_ok=True)
        stderr_log_output.write_text(f"miniprot not run: {reason}\n", encoding="utf-8")

    write_tsv(alignment_hits_output, alignment_rows(hits, full_coverage, partial_coverage), PHASE4B_ALIGNMENT_HIT_COLUMNS)
    write_tsv(
        summary_output,
        summary_rows(
            phase4_summary,
            queries,
            hits,
            ready,
            full_coverage,
            partial_coverage,
            eligible_genes=eligible_genes,
            genome_scope=genome_scope,
        ),
        PHASE4B_SUMMARY_COLUMNS,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Phase 4b genome-level validation for rescue candidates.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--phase4-summary", type=Path, required=True)
    parser.add_argument("--query-inventory", type=Path, required=True)
    parser.add_argument("--rescue-candidate-fasta", type=Path, required=True)
    parser.add_argument("--query-fasta-output", type=Path, required=True)
    parser.add_argument("--target-fasta-output", type=Path, required=True)
    parser.add_argument("--target-regions-output", type=Path, required=True)
    parser.add_argument("--preflight-output", type=Path, required=True)
    parser.add_argument("--raw-gff-output", type=Path, required=True)
    parser.add_argument("--stderr-log-output", type=Path, required=True)
    parser.add_argument("--alignment-hits-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    args = parser.parse_args()

    configure_logging()
    run_phase4b(
        args.config,
        args.phase4_summary,
        args.query_inventory,
        args.rescue_candidate_fasta,
        args.query_fasta_output,
        args.target_fasta_output,
        args.target_regions_output,
        args.preflight_output,
        args.raw_gff_output,
        args.stderr_log_output,
        args.alignment_hits_output,
        args.summary_output,
    )


if __name__ == "__main__":
    main()
