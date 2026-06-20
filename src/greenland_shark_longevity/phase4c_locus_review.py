"""Phase 4c locus review and TP53 targeted genome search."""

from __future__ import annotations

import argparse
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path

from .genome_validation import (
    QueryRecord,
    alignment_rows,
    extract_target_scaffolds,
    open_text,
    parse_fasta,
    parse_int,
    parse_miniprot_gff,
    resolve_miniprot,
    run_miniprot,
    write_query_fasta,
)
from .schemas import (
    DOMAIN_CHECK_TARGET_COLUMNS,
    DOMAIN_INTEGRITY_COLUMNS,
    PHASE4_RECIPROCAL_HIT_COLUMNS,
    PHASE4B_ALIGNMENT_HIT_COLUMNS,
    PHASE4B_QUERY_COLUMNS,
    PHASE4B_TARGET_REGION_COLUMNS,
    PHASE4C_GENE_REVIEW_COLUMNS,
    PHASE4C_LOCUS_REVIEW_COLUMNS,
    PHASE4C_TP53_CHUNK_PLAN_COLUMNS,
    PHASE4C_TP53_SUMMARY_COLUMNS,
    PHASE4C_TP53_TARGET_REGION_COLUMNS,
)
from .utils import configure_logging, read_tsv, read_yaml, write_tsv


NOT_ASSESSED = "NOT_ASSESSED"


@dataclass(frozen=True)
class PafTags:
    query_id: str
    target_seqid: str
    frameshifts: int
    stops: int
    donor_acceptor_changes: int
    donor_acceptor_disruptions: int
    cigar: str


@dataclass(frozen=True)
class FastaRecord:
    header: str
    seqid: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


def phase4c_config(config_path: Path) -> dict:
    config = read_yaml(config_path)
    return config.get("phase4c_locus_review", {})


def build_rescue_domain_targets(query_table: Path, query_fasta: Path) -> list[dict[str, str]]:
    rows = read_tsv(query_table, PHASE4B_QUERY_COLUMNS)
    targets: list[dict[str, str]] = []
    for row in rows:
        ready = row["query_status"] == "READY_FOR_GENOME_ALIGNMENT"
        targets.append(
            {
                "mechanism": row["mechanism"],
                "gene_symbol": row["gene_symbol"],
                "resource_id": row["source_resource_id"],
                "gene_id": row["original_protein_id"],
                "orthogroup_id": NOT_ASSESSED,
                "representative_protein_id": row["original_protein_id"],
                "representative_fasta_id": row["query_id"],
                "representative_length_aa": row["sequence_length_aa"],
                "candidate_fasta": str(query_fasta),
                "domain_check_status": "READY_FOR_DOMAIN_SCAN" if ready else "NOT_ASSESSED",
                "recommended_method": "HMMER/Pfam rescue-query domain screen",
                "required_validation": "Use domain support only as protein-model plausibility; inspect loci and orthology before biological interpretation.",
                "notes": "Phase 4c rescue-query domain target generated from Phase 4b miniprot query metadata.",
            }
        )
    return targets


def parse_paf_tags(raw_gff: Path) -> dict[str, PafTags]:
    tags_by_feature: dict[str, PafTags] = {}
    pending: PafTags | None = None
    if not raw_gff.exists():
        return tags_by_feature

    with raw_gff.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##PAF\t"):
                fields = line.split("\t")
                tag_values: dict[str, str] = {}
                for field in fields[13:]:
                    parts = field.split(":", 2)
                    if len(parts) == 3:
                        tag_values[parts[0]] = parts[2]
                pending = PafTags(
                    query_id=fields[1] if len(fields) > 1 else NOT_ASSESSED,
                    target_seqid=fields[6] if len(fields) > 6 else NOT_ASSESSED,
                    frameshifts=parse_int(tag_values.get("fs", "0")),
                    stops=parse_int(tag_values.get("st", "0")),
                    donor_acceptor_changes=parse_int(tag_values.get("da", "0")),
                    donor_acceptor_disruptions=parse_int(tag_values.get("do", "0")),
                    cigar=tag_values.get("cg", NOT_ASSESSED),
                )
                continue
            if pending is None or line.startswith("#") or not line:
                continue
            fields = line.split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "match", "protein_match"}:
                continue
            attributes = parse_attributes(fields[8])
            feature_id = attributes.get("ID")
            if feature_id:
                tags_by_feature[feature_id] = pending
            pending = None
    return tags_by_feature


def parse_attributes(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in value.split(";"):
        if not part:
            continue
        if "=" in part:
            key, attr_value = part.split("=", 1)
        else:
            key, attr_value = part, ""
        attrs[key] = attr_value
    return attrs


def overlap_fraction(a: dict[str, str], b: dict[str, str]) -> float:
    if a["target_seqid"] != b["target_seqid"]:
        return 0.0
    a_start = min(parse_int(a["target_start"]), parse_int(a["target_end"]))
    a_end = max(parse_int(a["target_start"]), parse_int(a["target_end"]))
    b_start = min(parse_int(b["target_start"]), parse_int(b["target_end"]))
    b_end = max(parse_int(b["target_start"]), parse_int(b["target_end"]))
    left = max(a_start, b_start)
    right = min(a_end, b_end)
    if right < left:
        return 0.0
    overlap = right - left + 1
    return overlap / max(1, min(a_end - a_start + 1, b_end - b_start + 1))


def cluster_hits(hits: list[dict[str, str]], min_overlap_fraction: float) -> list[list[dict[str, str]]]:
    clusters: list[list[dict[str, str]]] = []
    for hit in sorted(hits, key=lambda row: (row["gene_symbol"], row["target_seqid"], parse_int(row["target_start"]), parse_int(row["target_end"]), row["query_id"])):
        placed = False
        for cluster in clusters:
            if any(overlap_fraction(hit, other) >= min_overlap_fraction for other in cluster):
                cluster.append(hit)
                placed = True
                break
        if not placed:
            clusters.append([hit])
    return clusters


def scaffold_alias_map(target_regions_path: Path) -> dict[str, str]:
    aliases: dict[str, str] = {}
    if not target_regions_path.exists():
        return aliases
    for row in read_tsv(target_regions_path, PHASE4B_TARGET_REGION_COLUMNS):
        header = row["target_fasta_header"]
        if header and header != "NOT_FOUND":
            target_seqid = header.split()[0]
            aliases.setdefault(target_seqid, set()).add(row["requested_seqid"])  # type: ignore[union-attr]
    return {seqid: ",".join(sorted(values)) for seqid, values in aliases.items()}  # type: ignore[arg-type]


def domain_map(domain_integrity_path: Path) -> dict[str, dict[str, str]]:
    if not domain_integrity_path.exists():
        return {}
    return {
        row["representative_fasta_id"]: row
        for row in read_tsv(domain_integrity_path, DOMAIN_INTEGRITY_COLUMNS)
    }


def status_rank(status: str) -> int:
    ranks = {
        "DOMAIN_SUPPORTED": 4,
        "PARTIAL_DOMAIN": 3,
        "NO_EXPECTED_DOMAIN_DETECTED": 2,
        "NOT_ASSESSED": 1,
    }
    return ranks.get(status, 0)


def aggregate_domain_status(query_ids: list[str], domains: dict[str, dict[str, str]]) -> tuple[str, int, str]:
    rows = [domains[query_id] for query_id in query_ids if query_id in domains]
    if not rows:
        return (NOT_ASSESSED, 0, NOT_ASSESSED)
    best_status = max((row["domain_validation_status"] for row in rows), key=status_rank)
    supported = sum(1 for row in rows if row["domain_validation_status"] == "DOMAIN_SUPPORTED")
    labels = [
        f"{row['representative_fasta_id']}:{row['domain_validation_status']}:{row['best_pfam_accession']}:{row['best_pfam_name']}"
        for row in rows
    ]
    return (best_status, supported, ";".join(labels) if labels else NOT_ASSESSED)


def format_range(values: list[int]) -> str:
    if not values:
        return NOT_ASSESSED
    return f"{min(values)}-{max(values)}" if min(values) != max(values) else str(values[0])


def max_float(values: list[str]) -> str:
    parsed: list[float] = []
    for value in values:
        try:
            parsed.append(float(value))
        except ValueError:
            continue
    return f"{max(parsed):.6g}" if parsed else NOT_ASSESSED


def review_loci(
    alignment_hits_path: Path,
    raw_gff_path: Path,
    domain_integrity_path: Path,
    target_regions_path: Path,
    config_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    params = phase4c_config(config_path)
    genes = set(params.get("review_genes", ["H1F0", "FTH1B", "RAD51"]))
    min_overlap = float(params.get("min_overlap_fraction_for_same_locus", 0.50))
    high_status = {"HIGH_COVERAGE_NO_DISRUPTION", "POSSIBLE_DISRUPTION"}
    hits = [
        row
        for row in read_tsv(alignment_hits_path, PHASE4B_ALIGNMENT_HIT_COLUMNS)
        if row["gene_symbol"] in genes and row["alignment_status"] in high_status
    ]
    paf = parse_paf_tags(raw_gff_path)
    domains = domain_map(domain_integrity_path)
    mechanism_by_gene = {
        row["gene_symbol"]: row["mechanism"]
        for row in domains.values()
        if row.get("gene_symbol") and row.get("mechanism")
    }
    aliases = scaffold_alias_map(target_regions_path)

    locus_rows: list[dict[str, str]] = []
    gene_rows: list[dict[str, str]] = []
    for gene in sorted(genes):
        gene_hits = [row for row in hits if row["gene_symbol"] == gene]
        clusters = cluster_hits(gene_hits, min_overlap)
        for idx, cluster in enumerate(clusters, start=1):
            starts = [min(parse_int(row["target_start"]), parse_int(row["target_end"])) for row in cluster]
            ends = [max(parse_int(row["target_start"]), parse_int(row["target_end"])) for row in cluster]
            query_ids = sorted({row["query_id"] for row in cluster})
            original_ids = sorted({row["original_protein_id"] for row in cluster})
            feature_ids = [row["feature_id"] for row in cluster]
            paf_rows = [paf[feature_id] for feature_id in feature_ids if feature_id in paf]
            frameshifts = sum(item.frameshifts for item in paf_rows)
            stops = sum(item.stops for item in paf_rows)
            disruptions = frameshifts + stops
            domain_status, domain_supported_count, best_domain_hits = aggregate_domain_status(query_ids, domains)
            target_seqid = cluster[0]["target_seqid"]
            high_count = sum(1 for row in cluster if row["alignment_status"] == "HIGH_COVERAGE_NO_DISRUPTION")
            if disruptions:
                locus_status = "POSSIBLE_DISRUPTION_REQUIRES_MANUAL_REVIEW"
                artifact_risk = "high"
            elif domain_status == "DOMAIN_SUPPORTED" and high_count:
                locus_status = "HIGH_COVERAGE_DOMAIN_SUPPORTED_CANDIDATE_LOCUS"
                artifact_risk = "moderate"
            elif high_count:
                locus_status = "HIGH_COVERAGE_ALIGNMENT_DOMAIN_INCOMPLETE_OR_NOT_ASSESSED"
                artifact_risk = "moderate"
            else:
                locus_status = "ALIGNMENT_REQUIRES_MANUAL_REVIEW"
                artifact_risk = "high"

            locus_rows.append(
                {
                    "gene_symbol": gene,
                    "mechanism": mechanism_by_gene.get(gene, NOT_ASSESSED),
                    "locus_cluster_id": f"{gene}_L{idx:03d}",
                    "target_seqid": target_seqid,
                    "requested_scaffold_aliases": aliases.get(target_seqid, NOT_ASSESSED),
                    "cluster_start": str(min(starts)),
                    "cluster_end": str(max(ends)),
                    "strand_values": ",".join(sorted({row["strand"] for row in cluster})),
                    "supporting_query_count": str(len(query_ids)),
                    "supporting_query_ids": ";".join(query_ids),
                    "supporting_original_protein_ids": ";".join(original_ids),
                    "best_query_coverage": max_float([row["query_coverage"] for row in cluster]),
                    "best_identity": max_float([row["identity"] for row in cluster]),
                    "hit_count": str(len(cluster)),
                    "high_coverage_hit_count": str(high_count),
                    "cds_feature_count_range": format_range([parse_int(row["cds_feature_count"]) for row in cluster]),
                    "exon_feature_count_range": format_range([parse_int(row["exon_feature_count"]) for row in cluster]),
                    "paf_frameshift_count": str(frameshifts),
                    "paf_stop_count": str(stops),
                    "miniprot_disruption_status": "POTENTIAL_DISRUPTION_TAG" if disruptions else "NO_PARSED_FS_OR_ST_TAG",
                    "domain_validation_status": domain_status,
                    "domain_supported_query_count": str(domain_supported_count),
                    "best_domain_hits": best_domain_hits,
                    "overlap_cluster_size": str(len(cluster)),
                    "overlap_note": "Clustered because miniprot intervals overlap above the configured threshold." if len(cluster) > 1 else "Single non-overlapping candidate interval.",
                    "candidate_locus_status": locus_status,
                    "artifact_risk": artifact_risk,
                    "claim_language_guardrail": "Do not report validated duplication, function, adaptation, activation, inactivation, or absence from this review table alone.",
                    "required_validation": "Manual raw GFF/PAF review, exon/CDS inspection, domain review, synteny/local context, and cross-resource support.",
                    "notes": "Phase 4c locus review from miniprot candidate-scaffold output.",
                }
            )

        supported_loci = [row for row in locus_rows if row["gene_symbol"] == gene and row["domain_validation_status"] == "DOMAIN_SUPPORTED"]
        disrupted_loci = [row for row in locus_rows if row["gene_symbol"] == gene and row["miniprot_disruption_status"] == "POTENTIAL_DISRUPTION_TAG"]
        gene_loci = [row for row in locus_rows if row["gene_symbol"] == gene]
        gene_rows.append(
            {
                "gene_symbol": gene,
                "mechanism": mechanism_by_gene.get(gene, NOT_ASSESSED),
                "reviewed_locus_count": str(len(gene_loci)),
                "high_coverage_locus_count": str(sum(1 for row in gene_loci if parse_int(row["high_coverage_hit_count"]) > 0)),
                "domain_supported_locus_count": str(len(supported_loci)),
                "loci_with_disruption_flags": str(len(disrupted_loci)),
                "best_query_coverage": max_float([row["best_query_coverage"] for row in gene_loci]),
                "candidate_loci": ";".join(f"{row['target_seqid']}:{row['cluster_start']}-{row['cluster_end']}" for row in gene_loci) if gene_loci else NOT_ASSESSED,
                "phase4c_status": "MULTIPLE_CANDIDATE_LOCI_REQUIRE_MANUAL_REVIEW" if len(gene_loci) > 1 else ("ONE_CANDIDATE_LOCUS_REQUIRES_MANUAL_REVIEW" if gene_loci else "NO_REVIEWABLE_LOCUS"),
                "artifact_risk": "moderate" if gene_loci and not disrupted_loci else "high",
                "claim_language_guardrail": "Candidate loci are not validated duplications or biological mechanism evidence.",
                "required_validation": "Inspect raw alignments, domains, isoforms, local context, repeats if available, and cross-resource evidence before Phase 8 scoring.",
                "notes": "Phase 4c summarizes candidate-scaffold miniprot loci; this is not full-genome/cross-resource validation.",
            }
        )
    return locus_rows, gene_rows


def select_query(query_table: Path, query_fasta: Path, gene_symbol: str) -> QueryRecord | None:
    records = parse_fasta(query_fasta)
    for row in read_tsv(query_table, PHASE4B_QUERY_COLUMNS):
        if row["gene_symbol"] != gene_symbol or row["query_status"] != "READY_FOR_GENOME_ALIGNMENT":
            continue
        record = records.get(row["query_id"])
        if record is None:
            continue
        return QueryRecord(
            gene_symbol=row["gene_symbol"],
            mechanism=row["mechanism"],
            query_id=row["query_id"],
            query_source_type=row["query_source_type"],
            original_protein_id=row["original_protein_id"],
            original_description=row["original_description"],
            source_species_id=row["source_species_id"],
            source_resource_id=row["source_resource_id"],
            source_fasta=row["source_fasta"],
            sequence=record.sequence,
            selection_reason="Phase 4c targeted TP53 genome search query.",
        )
    return None


def iter_fasta_records(path: Path):
    header: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield FastaRecord(header=header, seqid=header.split()[0], sequence="".join(chunks))
                header = line[1:].strip()
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        yield FastaRecord(header=header, seqid=header.split()[0], sequence="".join(chunks))


def append_file(source: Path, target, header: str) -> None:
    target.write(header + "\n")
    if source.exists():
        with source.open("r", encoding="utf-8", errors="replace") as handle:
            shutil.copyfileobj(handle, target)


def write_chunk(path: Path, records: list[FastaRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            seq = record.sequence
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx : idx + 80] + "\n")


def run_tp53_chunked_search(
    config_path: Path,
    query_table: Path,
    query_fasta: Path,
    query_output: Path,
    chunk_plan_output: Path,
    raw_gff_output: Path,
    stderr_log_output: Path,
    alignment_hits_output: Path,
    summary_output: Path,
) -> None:
    config = read_yaml(config_path)
    phase4b = config.get("phase4b_genome_validation", {})
    params = phase4c_config(config_path)
    genome_fasta = Path(str(params.get("tp53_genome_fasta", phase4b.get("genome_fasta", ""))))
    miniprot = str(params.get("miniprot_executable", phase4b.get("miniprot_executable", "miniprot")))
    resolved_miniprot = resolve_miniprot(miniprot)
    threads = int(params.get("threads", phase4b.get("threads", 1)))
    max_chunk_bases = int(params.get("tp53_max_chunk_bases", 250_000_000))
    work_chunk = Path(str(params.get("tp53_work_chunk_fasta", "data/interim/annotation_rescue/tp53_active_chunk.fna")))
    chunk_gff = Path(str(params.get("tp53_work_chunk_gff", "data/interim/annotation_rescue/tp53_active_chunk.gff")))
    chunk_stderr = Path(str(params.get("tp53_work_chunk_stderr", "data/interim/annotation_rescue/tp53_active_chunk.stderr.log")))
    full_coverage = float(phase4b.get("min_query_coverage_for_intact_candidate", 0.70))
    partial_coverage = float(phase4b.get("min_query_coverage_for_partial_candidate", 0.20))

    query = select_query(query_table, query_fasta, "TP53")
    raw_gff_output.parent.mkdir(parents=True, exist_ok=True)
    stderr_log_output.parent.mkdir(parents=True, exist_ok=True)
    query_output.parent.mkdir(parents=True, exist_ok=True)
    plan_rows: list[dict[str, str]] = []
    hits = []

    if query is None or not genome_fasta.exists() or not resolved_miniprot:
        reason = "missing TP53 query, genome FASTA, or miniprot executable"
        if query:
            write_query_fasta(query_output, [query])
        else:
            query_output.write_text("", encoding="utf-8")
        raw_gff_output.write_text(f"##gff-version 3\n# tp53_search_not_run={reason}\n", encoding="utf-8")
        stderr_log_output.write_text(f"TP53 miniprot search not run: {reason}\n", encoding="utf-8")
        write_tsv(chunk_plan_output, plan_rows, PHASE4C_TP53_CHUNK_PLAN_COLUMNS)
        write_tsv(alignment_hits_output, [], PHASE4B_ALIGNMENT_HIT_COLUMNS)
        write_tsv(
            summary_output,
            [tp53_summary_row(query, [], plan_rows, ready=False, full_coverage=full_coverage, partial_coverage=partial_coverage)],
            PHASE4C_TP53_SUMMARY_COLUMNS,
        )
        return

    write_query_fasta(query_output, [query])
    query_lookup = {query.query_id: query}
    raw_gff_output.write_text("##gff-version 3\n# combined_phase4c_tp53_chunked_miniprot_output\n", encoding="utf-8")
    stderr_log_output.write_text("# combined_phase4c_tp53_chunked_miniprot_stderr\n", encoding="utf-8")

    chunk_records: list[FastaRecord] = []
    chunk_bases = 0
    chunk_id = 0

    def flush_chunk() -> None:
        nonlocal chunk_records, chunk_bases, chunk_id, hits
        if not chunk_records:
            return
        chunk_id += 1
        chunk_name = f"chunk_{chunk_id:04d}"
        write_chunk(work_chunk, chunk_records)
        logging.info("Running TP53 miniprot chunk %s with %s bases", chunk_name, chunk_bases)
        run_miniprot(resolved_miniprot, work_chunk, query_output, chunk_gff, chunk_stderr, threads, [])
        with raw_gff_output.open("a", encoding="utf-8") as raw_out:
            append_file(chunk_gff, raw_out, f"# phase4c_tp53_chunk={chunk_name}")
        with stderr_log_output.open("a", encoding="utf-8") as err_out:
            append_file(chunk_stderr, err_out, f"# phase4c_tp53_chunk={chunk_name}")
        chunk_hits = parse_miniprot_gff(chunk_gff, query_lookup)
        hits.extend(chunk_hits)
        plan_rows.append(
            {
                "chunk_id": chunk_name,
                "chunk_fasta": str(work_chunk),
                "sequence_count": str(len(chunk_records)),
                "total_bases": str(chunk_bases),
                "target_seqids": ";".join(record.seqid for record in chunk_records),
                "run_status": "MINIPROT_COMPLETED",
                "stderr_log": str(chunk_stderr),
                "raw_gff": str(chunk_gff),
                "notes": "Temporary chunk FASTA is overwritten between chunks; genome FASTA plus this plan define the search scope.",
            }
        )
        chunk_records = []
        chunk_bases = 0

    for record in iter_fasta_records(genome_fasta):
        if chunk_records and chunk_bases + record.length > max_chunk_bases:
            flush_chunk()
        chunk_records.append(record)
        chunk_bases += record.length
    flush_chunk()

    write_tsv(chunk_plan_output, plan_rows, PHASE4C_TP53_CHUNK_PLAN_COLUMNS)
    write_tsv(alignment_hits_output, alignment_rows(hits, full_coverage, partial_coverage), PHASE4B_ALIGNMENT_HIT_COLUMNS)
    write_tsv(
        summary_output,
        [tp53_summary_row(query, alignment_rows(hits, full_coverage, partial_coverage), plan_rows, ready=True, full_coverage=full_coverage, partial_coverage=partial_coverage)],
        PHASE4C_TP53_SUMMARY_COLUMNS,
    )


def tp53_summary_row(
    query: QueryRecord | None,
    hit_rows: list[dict[str, str]],
    plan_rows: list[dict[str, str]],
    ready: bool,
    full_coverage: float,
    partial_coverage: float,
) -> dict[str, str]:
    high = [row for row in hit_rows if float(row.get("query_coverage", "0") or 0) >= full_coverage and row.get("frameshift_or_stop_flag") != "True"]
    disrupted = [row for row in hit_rows if row.get("frameshift_or_stop_flag") == "True"]
    partial = [row for row in hit_rows if float(row.get("query_coverage", "0") or 0) >= partial_coverage]
    best_coverage = max([float(row.get("query_coverage", "0") or 0) for row in hit_rows], default=0.0)
    loci = [
        f"{row['target_seqid']}:{min(parse_int(row['target_start']), parse_int(row['target_end']))}-{max(parse_int(row['target_start']), parse_int(row['target_end']))}"
        for row in (high if high else partial)
    ]
    if not ready:
        status = "NOT_RUN_REQUIRES_INPUTS"
        classification = "not_assessable"
        notes = "TP53 targeted genome search did not run because required inputs were unavailable."
    elif high:
        status = "HIGH_COVERAGE_TP53_GENOME_ALIGNMENT"
        classification = "candidate_locus_requires_review"
        notes = "High-coverage TP53 protein-to-genome alignment detected; manual review is required."
    elif disrupted and best_coverage >= full_coverage:
        status = "POSSIBLE_DISRUPTED_TP53_ALIGNMENT"
        classification = "possible_disrupted_or_divergent_candidate_locus"
        notes = "A high-coverage TP53 protein-to-genome alignment was detected with parsed miniprot disruption flags; do not infer inactivation without manual review."
    elif partial:
        status = "PARTIAL_TP53_GENOME_ALIGNMENT"
        classification = "possible_fragment_or_divergent_locus"
        notes = "Only partial TP53 protein-to-genome alignment support was detected."
    else:
        status = "NO_TP53_ALIGNMENT_UNDER_CURRENT_QUERY"
        classification = "annotation_uncertainty"
        notes = "No qualifying TP53 alignment was detected with this comparator query; do not infer absence."
    return {
        "gene_symbol": "TP53",
        "query_id": query.query_id if query else NOT_ASSESSED,
        "query_source_type": query.query_source_type if query else NOT_ASSESSED,
        "query_length_aa": str(query.length) if query else NOT_ASSESSED,
        "search_scope": "chunked_full_genome",
        "chunks_scanned": str(len(plan_rows)),
        "total_bases_scanned": str(sum(parse_int(row["total_bases"]) for row in plan_rows)),
        "alignment_hit_count": str(len(hit_rows)),
        "high_coverage_hit_count": str(len(high)),
        "partial_hit_count": str(len(partial)),
        "best_query_coverage": f"{best_coverage:.6g}" if ready else NOT_ASSESSED,
        "candidate_loci": ";".join(loci) if loci else NOT_ASSESSED,
        "genome_search_status": status,
        "candidate_classification": classification,
        "artifact_risk": "high",
        "claim_language_guardrail": "Do not report TP53 absence, inactivation, activation, adaptation, or functional relevance from this search alone.",
        "required_validation": "Use additional TP53 queries, inspect raw miniprot alignments, compare resources, and validate domains/exons before interpretation.",
        "notes": notes,
    }


def parse_coordinate_seqid(value: str) -> str:
    parts = value.split(":")
    if len(parts) >= 3:
        return parts[1]
    return NOT_ASSESSED


def tp53_forward_target_scaffolds(reciprocal_hits_path: Path) -> tuple[dict[str, set[str]], list[dict[str, str]]]:
    scaffolds: set[str] = set()
    source_rows: list[dict[str, str]] = []
    for row in read_tsv(reciprocal_hits_path, PHASE4_RECIPROCAL_HIT_COLUMNS):
        if row["gene_symbol"] != "TP53":
            continue
        seqid = parse_coordinate_seqid(row["target_coordinate_summary"])
        if seqid == NOT_ASSESSED:
            continue
        scaffolds.add(seqid)
        source_rows.append(
            {
                "gene_symbol": "TP53",
                "source_table": str(reciprocal_hits_path),
                "source_protein_id": row["target_protein_id"],
                "source_description": row["target_description"],
                "source_coordinate": row["target_coordinate_summary"],
                "requested_seqid": seqid,
                "target_fasta_header": NOT_ASSESSED,
                "target_fasta_path": NOT_ASSESSED,
                "extraction_status": "PENDING",
                "notes": "TP53 targeted scaffold selected from Phase 4 forward-hit/reciprocal-audit coordinates; not reciprocal support.",
            }
        )
    return ({"TP53": scaffolds} if scaffolds else {}, source_rows)


def merge_tp53_target_rows(source_rows: list[dict[str, str]], extracted_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    extracted_by_seqid = {row["requested_seqid"]: row for row in extracted_rows}
    output: list[dict[str, str]] = []
    for row in source_rows:
        extracted = extracted_by_seqid.get(row["requested_seqid"], {})
        output.append(
            {
                **row,
                "target_fasta_header": extracted.get("target_fasta_header", NOT_ASSESSED),
                "target_fasta_path": extracted.get("target_fasta_path", NOT_ASSESSED),
                "extraction_status": extracted.get("extraction_status", "NOT_ASSESSED"),
                "notes": f"{row['notes']} Extraction status: {extracted.get('extraction_status', 'NOT_ASSESSED')}.",
            }
        )
    return output


def run_tp53_targeted_forward_search(
    config_path: Path,
    reciprocal_hits_path: Path,
    query_table: Path,
    query_fasta: Path,
    query_output: Path,
    target_fasta_output: Path,
    target_regions_output: Path,
    raw_gff_output: Path,
    stderr_log_output: Path,
    alignment_hits_output: Path,
    summary_output: Path,
) -> None:
    config = read_yaml(config_path)
    phase4b = config.get("phase4b_genome_validation", {})
    params = phase4c_config(config_path)
    genome_fasta = Path(str(params.get("tp53_genome_fasta", phase4b.get("genome_fasta", ""))))
    miniprot = str(params.get("miniprot_executable", phase4b.get("miniprot_executable", "miniprot")))
    resolved_miniprot = resolve_miniprot(miniprot)
    threads = int(params.get("threads", phase4b.get("threads", 1)))
    full_coverage = float(phase4b.get("min_query_coverage_for_intact_candidate", 0.70))
    partial_coverage = float(phase4b.get("min_query_coverage_for_partial_candidate", 0.20))
    query = select_query(query_table, query_fasta, "TP53")
    scaffold_map, source_rows = tp53_forward_target_scaffolds(reciprocal_hits_path)

    if query:
        write_query_fasta(query_output, [query])
    else:
        query_output.parent.mkdir(parents=True, exist_ok=True)
        query_output.write_text("", encoding="utf-8")

    if genome_fasta.exists() and scaffold_map:
        _, extracted_rows = extract_target_scaffolds(genome_fasta, target_fasta_output, scaffold_map)
    else:
        extracted_rows = []
        target_fasta_output.parent.mkdir(parents=True, exist_ok=True)
        target_fasta_output.write_text("", encoding="utf-8")
    write_tsv(target_regions_output, merge_tp53_target_rows(source_rows, extracted_rows), PHASE4C_TP53_TARGET_REGION_COLUMNS)

    ready = bool(query and resolved_miniprot and target_fasta_output.exists() and target_fasta_output.stat().st_size > 0)
    hits = []
    if ready:
        run_miniprot(resolved_miniprot, target_fasta_output, query_output, raw_gff_output, stderr_log_output, threads, [])
        hits = parse_miniprot_gff(raw_gff_output, {query.query_id: query})
    else:
        raw_gff_output.parent.mkdir(parents=True, exist_ok=True)
        stderr_log_output.parent.mkdir(parents=True, exist_ok=True)
        raw_gff_output.write_text("##gff-version 3\n# tp53_targeted_search_not_run=missing query, target FASTA, or miniprot\n", encoding="utf-8")
        stderr_log_output.write_text("TP53 targeted miniprot search not run: missing query, target FASTA, or miniprot\n", encoding="utf-8")

    rows = alignment_rows(hits, full_coverage, partial_coverage)
    write_tsv(alignment_hits_output, rows, PHASE4B_ALIGNMENT_HIT_COLUMNS)
    summary = tp53_summary_row(query, rows, [], ready=ready, full_coverage=full_coverage, partial_coverage=partial_coverage)
    summary["search_scope"] = "phase4_forward_hit_target_scaffolds"
    summary["chunks_scanned"] = "1" if ready else "0"
    summary["total_bases_scanned"] = str(sum(record.length for record in iter_fasta_records(target_fasta_output))) if target_fasta_output.exists() else "0"
    if ready and not rows:
        summary["notes"] = "No TP53 miniprot alignment was detected on Phase 4 forward-hit target scaffolds; do not infer absence."
    write_tsv(summary_output, [summary], PHASE4C_TP53_SUMMARY_COLUMNS)


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 4c locus review and TP53 targeted genome search.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    targets_parser = subparsers.add_parser("domain-targets")
    targets_parser.add_argument("--query-table", type=Path, required=True)
    targets_parser.add_argument("--query-fasta", type=Path, required=True)
    targets_parser.add_argument("--output", type=Path, required=True)

    review_parser = subparsers.add_parser("review-loci")
    review_parser.add_argument("--config", type=Path, required=True)
    review_parser.add_argument("--alignment-hits", type=Path, required=True)
    review_parser.add_argument("--raw-gff", type=Path, required=True)
    review_parser.add_argument("--domain-integrity", type=Path, required=True)
    review_parser.add_argument("--target-regions", type=Path, required=True)
    review_parser.add_argument("--locus-output", type=Path, required=True)
    review_parser.add_argument("--summary-output", type=Path, required=True)

    tp53_parser = subparsers.add_parser("tp53-search")
    tp53_parser.add_argument("--config", type=Path, required=True)
    tp53_parser.add_argument("--query-table", type=Path, required=True)
    tp53_parser.add_argument("--query-fasta", type=Path, required=True)
    tp53_parser.add_argument("--query-output", type=Path, required=True)
    tp53_parser.add_argument("--chunk-plan-output", type=Path, required=True)
    tp53_parser.add_argument("--raw-gff-output", type=Path, required=True)
    tp53_parser.add_argument("--stderr-log-output", type=Path, required=True)
    tp53_parser.add_argument("--alignment-hits-output", type=Path, required=True)
    tp53_parser.add_argument("--summary-output", type=Path, required=True)

    tp53_targeted_parser = subparsers.add_parser("tp53-targeted-forward-search")
    tp53_targeted_parser.add_argument("--config", type=Path, required=True)
    tp53_targeted_parser.add_argument("--reciprocal-hits", type=Path, required=True)
    tp53_targeted_parser.add_argument("--query-table", type=Path, required=True)
    tp53_targeted_parser.add_argument("--query-fasta", type=Path, required=True)
    tp53_targeted_parser.add_argument("--query-output", type=Path, required=True)
    tp53_targeted_parser.add_argument("--target-fasta-output", type=Path, required=True)
    tp53_targeted_parser.add_argument("--target-regions-output", type=Path, required=True)
    tp53_targeted_parser.add_argument("--raw-gff-output", type=Path, required=True)
    tp53_targeted_parser.add_argument("--stderr-log-output", type=Path, required=True)
    tp53_targeted_parser.add_argument("--alignment-hits-output", type=Path, required=True)
    tp53_targeted_parser.add_argument("--summary-output", type=Path, required=True)

    args = parser.parse_args()
    configure_logging()
    if args.command == "domain-targets":
        write_tsv(args.output, build_rescue_domain_targets(args.query_table, args.query_fasta), DOMAIN_CHECK_TARGET_COLUMNS)
    elif args.command == "review-loci":
        locus_rows, summary_rows = review_loci(
            args.alignment_hits,
            args.raw_gff,
            args.domain_integrity,
            args.target_regions,
            args.config,
        )
        write_tsv(args.locus_output, locus_rows, PHASE4C_LOCUS_REVIEW_COLUMNS)
        write_tsv(args.summary_output, summary_rows, PHASE4C_GENE_REVIEW_COLUMNS)
    elif args.command == "tp53-search":
        run_tp53_chunked_search(
            args.config,
            args.query_table,
            args.query_fasta,
            args.query_output,
            args.chunk_plan_output,
            args.raw_gff_output,
            args.stderr_log_output,
            args.alignment_hits_output,
            args.summary_output,
        )
    elif args.command == "tp53-targeted-forward-search":
        run_tp53_targeted_forward_search(
            args.config,
            args.reciprocal_hits,
            args.query_table,
            args.query_fasta,
            args.query_output,
            args.target_fasta_output,
            args.target_regions_output,
            args.raw_gff_output,
            args.stderr_log_output,
            args.alignment_hits_output,
            args.summary_output,
        )
    else:
        raise ValueError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
