"""Phase 4 targeted annotation rescue for unresolved high-priority genes."""

from __future__ import annotations

import argparse
import logging
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from pyhmmer import easel
from pyhmmer.hmmer import phmmer

from .candidate_orthofinder import TARGET_RESOURCE_ID, TARGET_SPECIES_ID
from .schemas import (
    ORTHOFINDER_INPUT_MANIFEST_COLUMNS,
    PHASE4_FORWARD_HIT_COLUMNS,
    PHASE4_GENOME_RESCUE_PLAN_COLUMNS,
    PHASE4_RECIPROCAL_HIT_COLUMNS,
    PHASE4_RESCUE_QUERY_COLUMNS,
    PHASE4_RESCUE_SUMMARY_COLUMNS,
    REFERENCE_GENE_COORDINATE_COLUMNS,
    RESCUE_TARGET_COLUMNS,
)
from .utils import configure_logging, read_tsv, read_yaml, write_tsv


ANNOTATION_UNCERTAINTY = "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH"

DEFAULT_PATTERNS = {
    "H1F0": {
        "include": ["histone h1.0", "histone h1.0-like", "histone h1.0-b-like"],
        "exclude": ["h1.01", "h1.10", "h1.1", "h1.2", "h1.3", "h1.4", "h1.5"],
    },
    "FTH1B": {
        "include": ["ferritin heavy chain b", "ferritin heavy chain b-like"],
        "exclude": ["light chain", "middle subunit", "lower subunit", "oocyte"],
    },
    "TP53": {
        "include": ["cellular tumor antigen p53"],
        "exclude": ["binding protein", "inducible", "regulated", "tp53rk", "p53 response", "p53-target"],
    },
    "RAD51": {
        "include": ["dna repair protein rad51 homolog 1"],
        "exclude": ["homolog 2", "homolog 3", "homolog 4", "associated"],
    },
}


@dataclass(frozen=True)
class FastaRecord:
    protein_id: str
    description: str
    sequence: str
    source_fasta: str
    species_id: str
    resource_id: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def header(self) -> str:
        return f"{self.protein_id} {self.description}".strip()


@dataclass(frozen=True)
class SearchHit:
    query: FastaRecord
    target: FastaRecord
    evalue: float
    score: float
    domain_i_evalue: float
    domain_score: float
    query_coverage: float
    target_coverage: float
    env_from: int
    env_to: int
    passes_filters: bool


def phase4_config(config_path: Path) -> dict:
    config = read_yaml(config_path)
    params = config.get("phase4_annotation_rescue", {})
    if not isinstance(params, dict):
        raise ValueError("config/phase4_annotation_rescue must be a mapping")
    return params


def parse_fasta(path: Path, species_id: str, resource_id: str) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    current_header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append(record_from_header(current_header, chunks, path, species_id, resource_id))
                current_header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if current_header is not None:
        records.append(record_from_header(current_header, chunks, path, species_id, resource_id))
    return records


def record_from_header(header: str, chunks: list[str], path: Path, species_id: str, resource_id: str) -> FastaRecord:
    protein_id, _, description = header.partition(" ")
    sequence = sanitize_sequence("".join(chunks))
    return FastaRecord(
        protein_id=protein_id,
        description=description,
        sequence=sequence,
        source_fasta=str(path),
        species_id=species_id,
        resource_id=resource_id,
    )


def sanitize_sequence(sequence: str) -> str:
    cleaned = []
    for residue in sequence.upper():
        if residue == "*":
            continue
        cleaned.append(residue if "A" <= residue <= "Z" else "X")
    return "".join(cleaned)


def digital_sequence(record: FastaRecord, alphabet: easel.Alphabet) -> easel.DigitalSequence:
    return easel.TextSequence(
        name=record.protein_id,
        description=record.description,
        sequence=record.sequence,
    ).digitize(alphabet)


def manifest_records(manifest_path: Path) -> list[dict[str, str]]:
    return [row for row in read_tsv(manifest_path, ORTHOFINDER_INPUT_MANIFEST_COLUMNS) if row["orthofinder_ready"] == "True"]


def load_proteomes(manifest_path: Path) -> dict[str, list[FastaRecord]]:
    proteomes: dict[str, list[FastaRecord]] = {}
    for row in manifest_records(manifest_path):
        path = Path(row["orthofinder_input_path"])
        if not path.exists():
            raise FileNotFoundError(f"Staged OrthoFinder FASTA is missing: {path}")
        proteomes[row["species_id"]] = parse_fasta(path, row["species_id"], row["resource_id"])
    return proteomes


def patterns_for(gene_symbol: str, config_patterns: dict | None = None) -> dict[str, list[str]]:
    merged = DEFAULT_PATTERNS.get(gene_symbol.upper(), {"include": [gene_symbol.lower()], "exclude": []})
    if config_patterns and gene_symbol in config_patterns:
        override = config_patterns[gene_symbol]
        if isinstance(override, dict):
            merged = {
                "include": [str(value).lower() for value in override.get("include", merged["include"])],
                "exclude": [str(value).lower() for value in override.get("exclude", merged["exclude"])],
            }
    return merged


def record_matches(record: FastaRecord, gene_symbol: str, config_patterns: dict | None = None) -> bool:
    text = record.header.lower()
    patterns = patterns_for(gene_symbol, config_patterns)
    return any(pattern in text for pattern in patterns["include"]) and not any(
        pattern in text for pattern in patterns["exclude"]
    )


def quality_key(record: FastaRecord) -> tuple[int, int, str]:
    text = record.header.lower()
    low_quality = int("low quality" in text or "partial" in text)
    return (low_quality, -record.length, record.protein_id)


def unresolved_targets(rescue_targets_path: Path, high_priority_genes: list[str]) -> list[dict[str, str]]:
    rows = read_tsv(rescue_targets_path, RESCUE_TARGET_COLUMNS)
    allowed = {gene.upper() for gene in high_priority_genes}
    return [
        row
        for row in rows
        if row["gene_symbol"].upper() in allowed and row["mapping_status"] == ANNOTATION_UNCERTAINTY
    ]


def select_query_records(
    targets: list[dict[str, str]],
    proteomes: dict[str, list[FastaRecord]],
    max_per_species: int,
    config_patterns: dict | None = None,
) -> tuple[list[dict[str, str]], dict[str, list[FastaRecord]]]:
    query_rows: list[dict[str, str]] = []
    selected: dict[str, list[FastaRecord]] = {}
    for target in targets:
        gene_symbol = target["gene_symbol"].upper()
        selected[gene_symbol] = []
        found_any = False
        for species_id, records in sorted(proteomes.items()):
            if species_id == TARGET_SPECIES_ID:
                continue
            matches = sorted(
                [record for record in records if record_matches(record, gene_symbol, config_patterns)],
                key=quality_key,
            )
            for rank, record in enumerate(matches):
                status = "SELECTED_FOR_RESCUE_QUERY" if rank < max_per_species else "MATCHED_BUT_EXCLUDED_BY_QUERY_CAP"
                if rank < max_per_species:
                    selected[gene_symbol].append(record)
                    found_any = True
                query_rows.append(
                    {
                        "gene_symbol": gene_symbol,
                        "mechanism": target["mechanism"],
                        "query_species_id": record.species_id,
                        "query_resource_id": record.resource_id,
                        "query_protein_id": record.protein_id,
                        "query_description": record.description,
                        "query_length_aa": str(record.length),
                        "source_fasta": record.source_fasta,
                        "selection_status": status,
                        "selection_rule": "description pattern match with low-quality/partial records deprioritized",
                        "notes": "Selected only as a query for annotation rescue; query annotation is not biological evidence.",
                    }
                )
        if not found_any:
            query_rows.append(
                {
                    "gene_symbol": gene_symbol,
                    "mechanism": target["mechanism"],
                    "query_species_id": "NOT_ASSESSED",
                    "query_resource_id": "NOT_ASSESSED",
                    "query_protein_id": "NOT_ASSESSED",
                    "query_description": "NOT_ASSESSED",
                    "query_length_aa": "NOT_ASSESSED",
                    "source_fasta": "NOT_ASSESSED",
                    "selection_status": "NO_QUERY_FOUND",
                    "selection_rule": "description pattern match with low-quality/partial records deprioritized",
                    "notes": "No comparator query met the configured pattern filters.",
                }
            )
    return query_rows, selected


def best_domain_values(hit, query_length: int) -> tuple[float, float, float, float, int, int]:
    if len(hit.domains) == 0 or hit.best_domain is None:
        return hit.evalue, 0.0, 0.0, 0.0, 0, 0
    domain = hit.best_domain
    alignment = domain.alignment
    query_span = max(0, alignment.hmm_to - alignment.hmm_from + 1)
    target_span = max(0, alignment.target_to - alignment.target_from + 1)
    query_coverage = query_span / query_length if query_length else 0.0
    target_coverage = target_span / hit.length if hit.length else 0.0
    return (
        domain.i_evalue,
        domain.score,
        query_coverage,
        target_coverage,
        domain.env_from,
        domain.env_to,
    )


def run_search(
    queries: list[FastaRecord],
    database: list[FastaRecord],
    cpus: int,
    max_evalue: float,
    min_query_coverage: float,
    min_target_coverage: float,
    max_hits_per_query: int,
) -> list[SearchHit]:
    if not queries or not database:
        return []
    alphabet = easel.Alphabet.amino()
    digital_queries = [digital_sequence(record, alphabet) for record in queries]
    digital_database = [digital_sequence(record, alphabet) for record in database]
    database_by_id = {record.protein_id: record for record in database}
    hits: list[SearchHit] = []
    for query, top_hits in zip(queries, phmmer(digital_queries, digital_database, cpus=cpus, E=max_evalue, domE=max_evalue)):
        ranked_hits = list(top_hits)[:max_hits_per_query]
        for hit in ranked_hits:
            target_name = decode_name(hit.name)
            target_record = database_by_id.get(target_name)
            if target_record is None:
                continue
            domain_i_evalue, domain_score, query_coverage, target_coverage, env_from, env_to = best_domain_values(
                hit, query.length
            )
            passes = (
                hit.evalue <= max_evalue
                and domain_i_evalue <= max_evalue
                and query_coverage >= min_query_coverage
                and target_coverage >= min_target_coverage
            )
            hits.append(
                SearchHit(
                    query=query,
                    target=target_record,
                    evalue=hit.evalue,
                    score=hit.score,
                    domain_i_evalue=domain_i_evalue,
                    domain_score=domain_score,
                    query_coverage=query_coverage,
                    target_coverage=target_coverage,
                    env_from=env_from,
                    env_to=env_to,
                    passes_filters=passes,
                )
            )
    return hits


def decode_name(value) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def fmt_float(value: float) -> str:
    return f"{value:.6g}"


def coordinate_lookup(coordinates_path: Path) -> dict[str, str]:
    rows = read_tsv(coordinates_path, REFERENCE_GENE_COORDINATE_COLUMNS)
    lookup: dict[str, str] = {}
    for row in rows:
        if row["resource_id"] != TARGET_RESOURCE_ID or row["parse_status"] != "PARSED_FROM_GFF":
            continue
        lookup[row["gene_id"]] = f"{row['gene_id']}:{row['seqid']}:{row['start']}-{row['end']}:{row['strand']}"
    return lookup


def gene_id_from_protein(protein_id: str) -> str:
    match = re.search(r"(gs_\d+)-P\d+", protein_id)
    if match:
        return f"gene-{match.group(1)}"
    return "NOT_ASSESSED"


def forward_rows(hits: Iterable[SearchHit], max_evalue: float, min_query_coverage: float, min_target_coverage: float) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for hit in hits:
        rows.append(
            {
                "gene_symbol": infer_gene_from_query(hit.query),
                "query_species_id": hit.query.species_id,
                "query_resource_id": hit.query.resource_id,
                "query_protein_id": hit.query.protein_id,
                "query_description": hit.query.description,
                "query_length_aa": str(hit.query.length),
                "target_species_id": hit.target.species_id,
                "target_resource_id": hit.target.resource_id,
                "target_protein_id": hit.target.protein_id,
                "target_description": hit.target.description,
                "target_length_aa": str(hit.target.length),
                "full_sequence_evalue": fmt_float(hit.evalue),
                "full_sequence_bitscore": fmt_float(hit.score),
                "best_domain_i_evalue": fmt_float(hit.domain_i_evalue),
                "best_domain_bitscore": fmt_float(hit.domain_score),
                "query_coverage": fmt_float(hit.query_coverage),
                "target_coverage": fmt_float(hit.target_coverage),
                "target_env_from": str(hit.env_from),
                "target_env_to": str(hit.env_to),
                "passes_similarity_filters": str(hit.passes_filters),
                "filter_rule": f"evalue <= {max_evalue:g}; query_coverage >= {min_query_coverage:g}; target_coverage >= {min_target_coverage:g}",
                "notes": "Forward protein similarity hit; not a gene model validation by itself.",
            }
        )
    return rows


def infer_gene_from_query(record: FastaRecord) -> str:
    # Query records are tagged externally during search orchestration.
    match = re.match(r"([A-Z0-9]+)\|", record.description)
    return match.group(1) if match else "NOT_ASSESSED"


def tag_query(record: FastaRecord, gene_symbol: str) -> FastaRecord:
    return FastaRecord(
        protein_id=record.protein_id,
        description=f"{gene_symbol}|{record.description}",
        sequence=record.sequence,
        source_fasta=record.source_fasta,
        species_id=record.species_id,
        resource_id=record.resource_id,
    )


def reciprocal_rows(
    forward_hits: list[SearchHit],
    proteomes: dict[str, list[FastaRecord]],
    selected_queries: dict[str, list[FastaRecord]],
    coordinates: dict[str, str],
    cpus: int,
    max_evalue: float,
    min_query_coverage: float,
    min_target_coverage: float,
    config_patterns: dict | None = None,
) -> tuple[list[dict[str, str]], dict[str, list[SearchHit]]]:
    rows: list[dict[str, str]] = []
    supported: dict[str, list[SearchHit]] = {}
    seen: set[tuple[str, str, str]] = set()
    for hit in forward_hits:
        if not hit.passes_filters:
            continue
        gene_symbol = infer_gene_from_query(hit.query)
        key = (gene_symbol, hit.target.protein_id, hit.query.species_id)
        if key in seen:
            continue
        seen.add(key)
        source_db = proteomes.get(hit.query.species_id, [])
        reciprocal_hits = run_search(
            [hit.target],
            source_db,
            cpus=cpus,
            max_evalue=max_evalue,
            min_query_coverage=min_query_coverage,
            min_target_coverage=min_target_coverage,
            max_hits_per_query=1,
        )
        reciprocal = reciprocal_hits[0] if reciprocal_hits else None
        selected_query_ids = {query.protein_id for query in selected_queries.get(gene_symbol, [])}
        if reciprocal:
            pattern_match = record_matches(reciprocal.target, gene_symbol, config_patterns)
            selected_match = reciprocal.target.protein_id in selected_query_ids
            if selected_match:
                support = "RECIPROCAL_SELECTED_QUERY_TOP_HIT"
            elif pattern_match:
                support = "RECIPROCAL_GENE_PATTERN_TOP_HIT"
            else:
                support = "NOT_RECIPROCAL"
        else:
            pattern_match = False
            selected_match = False
            support = "NO_RECIPROCAL_HIT"
        if support != "NOT_RECIPROCAL" and support != "NO_RECIPROCAL_HIT":
            supported.setdefault(gene_symbol, []).append(hit)
        target_gene_id = gene_id_from_protein(hit.target.protein_id)
        rows.append(
            {
                "gene_symbol": gene_symbol,
                "target_protein_id": hit.target.protein_id,
                "target_description": hit.target.description,
                "target_gene_id": target_gene_id,
                "target_coordinate_summary": coordinates.get(target_gene_id, "NO_COORDINATES"),
                "source_species_id": hit.query.species_id,
                "source_resource_id": hit.query.resource_id,
                "query_protein_id": hit.query.protein_id,
                "reciprocal_top_protein_id": reciprocal.target.protein_id if reciprocal else "NOT_ASSESSED",
                "reciprocal_top_description": reciprocal.target.description if reciprocal else "NOT_ASSESSED",
                "reciprocal_top_evalue": fmt_float(reciprocal.evalue) if reciprocal else "NOT_ASSESSED",
                "reciprocal_top_bitscore": fmt_float(reciprocal.score) if reciprocal else "NOT_ASSESSED",
                "reciprocal_gene_pattern_match": str(pattern_match),
                "reciprocal_selected_query_match": str(selected_match),
                "reciprocal_support_status": support,
                "notes": "Reciprocal protein similarity support; still requires domain and genome-context review.",
            }
        )
    return rows, supported


def summary_rows(
    targets: list[dict[str, str]],
    selected_queries: dict[str, list[FastaRecord]],
    forward_hits: list[SearchHit],
    reciprocal_supported: dict[str, list[SearchHit]],
    coordinates: dict[str, str],
    genome_fasta: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for target in targets:
        gene_symbol = target["gene_symbol"].upper()
        passing_forward = [
            hit for hit in forward_hits if infer_gene_from_query(hit.query) == gene_symbol and hit.passes_filters
        ]
        supported_hits = reciprocal_supported.get(gene_symbol, [])
        candidate_gene_ids = sorted({gene_id_from_protein(hit.target.protein_id) for hit in supported_hits})
        candidate_proteins = sorted({hit.target.protein_id for hit in supported_hits})
        candidate_products = sorted({hit.target.description for hit in supported_hits})
        coordinate_summary = ";".join(coordinates.get(gene_id, f"{gene_id}:NO_COORDINATES") for gene_id in candidate_gene_ids)
        if supported_hits:
            protein_status = "PROTEIN_LEVEL_RESCUE_CANDIDATE"
            classification = "Tier 2 candidate rescue input"
            artifact_risk = "moderate_to_high"
            notes = "Reciprocal protein similarity supports one or more focal annotation candidates, but this is not final annotation validation."
        elif selected_queries.get(gene_symbol):
            protein_status = "NO_RECIPROCAL_PROTEIN_RESCUE"
            classification = "Artifact/uncertain"
            artifact_risk = "high"
            notes = "Queries were available, but no reciprocal protein-level rescue candidate passed the configured filters."
        else:
            protein_status = "QUERY_NOT_AVAILABLE"
            classification = "NOT_ASSESSED"
            artifact_risk = "not_assessable"
            notes = "No comparator query sequence met the configured selection rules."

        rows.append(
            {
                "gene_symbol": gene_symbol,
                "mechanism": target["mechanism"],
                "initial_mapping_status": target["mapping_status"],
                "query_count": str(len(selected_queries.get(gene_symbol, []))),
                "passing_forward_hit_count": str(len(passing_forward)),
                "reciprocal_supported_hit_count": str(len(supported_hits)),
                "candidate_gene_ids": ",".join(candidate_gene_ids) if candidate_gene_ids else "NOT_ASSESSED",
                "candidate_protein_ids": ",".join(candidate_proteins) if candidate_proteins else "NOT_ASSESSED",
                "candidate_products": " | ".join(candidate_products) if candidate_products else "NOT_ASSESSED",
                "candidate_coordinate_summary": coordinate_summary if coordinate_summary else "NOT_ASSESSED",
                "protein_rescue_status": protein_status,
                "genome_alignment_status": "READY_FOR_MINIPROT" if genome_fasta.exists() else "GENOME_FASTA_MISSING_NOT_RUN",
                "classification": classification,
                "artifact_risk": artifact_risk,
                "claim_language_guardrail": "Do not report absent, inactivated, activated, adaptive, causal, or validated duplication language from this rescue table alone.",
                "required_validation": "Run domain checks on rescued candidates, inspect separable loci/isoforms, and run protein-to-genome alignment when local genome FASTA is available.",
                "notes": notes,
            }
        )
    return rows


def genome_plan_rows(targets: list[dict[str, str]], selected_queries: dict[str, list[FastaRecord]], genome_fasta: Path) -> list[dict[str, str]]:
    miniprot_path = shutil.which("miniprot")
    rows: list[dict[str, str]] = []
    for target in targets:
        gene_symbol = target["gene_symbol"].upper()
        rows.append(
            {
                "gene_symbol": gene_symbol,
                "recommended_method": "miniprot protein-to-genome alignment after reciprocal protein evidence identifies query proteins.",
                "genome_fasta": str(genome_fasta),
                "genome_fasta_status": "LOCAL_AVAILABLE" if genome_fasta.exists() else "MISSING_LOCAL",
                "miniprot_status": "AVAILABLE_ON_PATH" if miniprot_path else "NOT_AVAILABLE_ON_PATH",
                "protein_query_status": "LOCAL_QUERY_AVAILABLE" if selected_queries.get(gene_symbol) else "NO_QUERY_AVAILABLE",
                "run_status": "NOT_RUN_REQUIRES_INPUTS" if not genome_fasta.exists() or not miniprot_path else "READY_NOT_RUN",
                "blockers": "Genome FASTA and miniprot are required; do not infer absence from missing local genome alignment.",
                "notes": "Protein-to-genome rescue is intentionally not run until required local inputs are present.",
            }
        )
    return rows


def write_fasta(path: Path, records: Iterable[FastaRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    with path.open("w", encoding="utf-8") as handle:
        for record in sorted(records, key=lambda item: item.protein_id):
            if record.protein_id in seen:
                continue
            seen.add(record.protein_id)
            handle.write(f">{record.protein_id} {record.description}\n")
            for idx in range(0, len(record.sequence), 80):
                handle.write(record.sequence[idx : idx + 80] + "\n")


def run_phase4(
    config_path: Path,
    rescue_targets_path: Path,
    manifest_path: Path,
    coordinates_path: Path,
    query_inventory_output: Path,
    forward_hits_output: Path,
    reciprocal_hits_output: Path,
    summary_output: Path,
    genome_plan_output: Path,
    candidate_fasta_output: Path,
) -> None:
    config = read_yaml(config_path)
    params = config.get("phase4_annotation_rescue", {})
    high_priority_genes = [str(gene).upper() for gene in params.get("high_priority_genes", ["H1F0", "FTH1B", "TP53", "RAD51"])]
    targets = unresolved_targets(rescue_targets_path, high_priority_genes)
    proteomes = load_proteomes(manifest_path)
    target_proteome = proteomes.get(TARGET_SPECIES_ID, [])
    if not target_proteome:
        raise ValueError(f"Target proteome for species_id={TARGET_SPECIES_ID} was not found in {manifest_path}")

    query_rows, selected_queries = select_query_records(
        targets,
        proteomes,
        max_per_species=int(params.get("max_queries_per_gene_per_species", 3)),
        config_patterns=params.get("query_patterns"),
    )
    tagged_queries = [
        tag_query(record, gene_symbol)
        for gene_symbol, records in selected_queries.items()
        for record in records
    ]
    cpus = int(params.get("threads", 1))
    max_evalue = float(params.get("max_evalue", 1e-5))
    min_query_coverage = float(params.get("min_query_coverage", 0.40))
    min_target_coverage = float(params.get("min_target_coverage", 0.40))
    max_hits = int(params.get("max_forward_hits_per_query", 5))
    logging.info("Running Phase 4 forward phmmer rescue with %s query proteins", len(tagged_queries))
    forward = run_search(
        tagged_queries,
        target_proteome,
        cpus=cpus,
        max_evalue=max_evalue,
        min_query_coverage=min_query_coverage,
        min_target_coverage=min_target_coverage,
        max_hits_per_query=max_hits,
    )
    coordinates = coordinate_lookup(coordinates_path)
    reciprocal, supported = reciprocal_rows(
        forward,
        proteomes,
        selected_queries,
        coordinates,
        cpus=cpus,
        max_evalue=max_evalue,
        min_query_coverage=min_query_coverage,
        min_target_coverage=min_target_coverage,
        config_patterns=params.get("query_patterns"),
    )
    genome_fasta = Path(str(config.get("reference_inputs", {}).get("assembly_fasta", "NOT_ASSESSED")))
    summary = summary_rows(targets, selected_queries, forward, supported, coordinates, genome_fasta)
    genome_plan = genome_plan_rows(targets, selected_queries, genome_fasta)
    candidate_records = [hit.target for hits in supported.values() for hit in hits]

    write_tsv(query_inventory_output, query_rows, PHASE4_RESCUE_QUERY_COLUMNS)
    write_tsv(
        forward_hits_output,
        forward_rows(forward, max_evalue, min_query_coverage, min_target_coverage),
        PHASE4_FORWARD_HIT_COLUMNS,
    )
    write_tsv(reciprocal_hits_output, reciprocal, PHASE4_RECIPROCAL_HIT_COLUMNS)
    write_tsv(summary_output, summary, PHASE4_RESCUE_SUMMARY_COLUMNS)
    write_tsv(genome_plan_output, genome_plan, PHASE4_GENOME_RESCUE_PLAN_COLUMNS)
    write_fasta(candidate_fasta_output, candidate_records)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run targeted Phase 4 annotation rescue on unresolved candidates.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--rescue-targets", type=Path, required=True)
    parser.add_argument("--orthofinder-input-manifest", type=Path, required=True)
    parser.add_argument("--gene-coordinates", type=Path, required=True)
    parser.add_argument("--query-inventory-output", type=Path, required=True)
    parser.add_argument("--forward-hits-output", type=Path, required=True)
    parser.add_argument("--reciprocal-hits-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--genome-plan-output", type=Path, required=True)
    parser.add_argument("--candidate-fasta-output", type=Path, required=True)
    args = parser.parse_args()

    configure_logging()
    run_phase4(
        args.config,
        args.rescue_targets,
        args.orthofinder_input_manifest,
        args.gene_coordinates,
        args.query_inventory_output,
        args.forward_hits_output,
        args.reciprocal_hits_output,
        args.summary_output,
        args.genome_plan_output,
        args.candidate_fasta_output,
    )


if __name__ == "__main__":
    main()
