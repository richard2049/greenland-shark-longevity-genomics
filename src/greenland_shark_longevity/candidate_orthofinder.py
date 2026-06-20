"""Map curated candidates onto real OrthoFinder orthogroups."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from .candidate_panels import iter_candidates
from .schemas import (
    COPY_NUMBER_COLUMNS,
    DUPLICATION_AUDIT_COLUMNS,
    ORTHOFINDER_INPUT_MANIFEST_COLUMNS,
    ORTHOGROUP_GENE_COUNT_LONG_COLUMNS,
    REFERENCE_GENE_COORDINATE_COLUMNS,
)
from .utils import configure_logging, join_values, read_tsv, read_yaml, write_tsv


TARGET_SPECIES_ID = "smic"
TARGET_RESOURCE_ID = "SMIC_TOKYO_GENOME_2025"


@dataclass(frozen=True)
class ProteinAnnotation:
    protein_id: str
    gene_id: str
    gene_symbol: str
    locus_tag: str
    product: str


def parse_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for part in value.split(";"):
        if not part or "=" not in part:
            continue
        key, raw = part.split("=", 1)
        attributes[key] = unquote(raw)
    return attributes


def normalise_label(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def candidate_aliases(candidate: dict[str, object]) -> set[str]:
    aliases = {normalise_label(str(candidate["gene_symbol"]))}
    for synonym in candidate.get("synonyms", []):
        synonym_text = str(synonym).strip()
        if not synonym_text:
            continue
        # Use compact aliases such as H1-0, CHK2, FTH1b, or p53. Descriptive
        # phrases are kept out of symbol matching to avoid product-name hits.
        if " " not in synonym_text and len(synonym_text) <= 20:
            aliases.add(normalise_label(synonym_text))
    return aliases


def orthofinder_identifier_variants(protein_id: str) -> set[str]:
    variants = {protein_id}
    variants.add(protein_id.replace(":", "_"))
    variants.add(protein_id.replace("|WGS:", "|WGS_"))
    variants.add(protein_id.replace("|", "_").replace(":", "_"))
    return variants


def _gene_id_from_locus_tag(locus_tag: str, fallback: str) -> str:
    if locus_tag:
        return f"gene-{locus_tag}"
    return fallback or "NOT_ASSESSED"


def parse_gff_protein_annotations(gff_path: Path) -> list[ProteinAnnotation]:
    annotations: dict[str, ProteinAnnotation] = {}
    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            attrs = parse_attributes(parts[8])
            protein_id = attrs.get("protein_id") or attrs.get("Name") or ""
            if not protein_id or protein_id in annotations:
                continue
            locus_tag = attrs.get("locus_tag", "")
            gene_symbol = attrs.get("gene") or attrs.get("Name") or locus_tag or "NOT_ASSESSED"
            annotations[protein_id] = ProteinAnnotation(
                protein_id=protein_id,
                gene_id=_gene_id_from_locus_tag(locus_tag, attrs.get("Parent", "")),
                gene_symbol=gene_symbol,
                locus_tag=locus_tag,
                product=attrs.get("product", "NOT_ASSESSED"),
            )
    return list(annotations.values())


def parse_orthogroups(orthogroups_path: Path) -> tuple[dict[str, dict[str, list[str]]], dict[str, tuple[str, str]]]:
    orthogroups: dict[str, dict[str, list[str]]] = {}
    protein_to_group: dict[str, tuple[str, str]] = {}
    with orthogroups_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "Orthogroup" not in reader.fieldnames:
            raise ValueError(f"{orthogroups_path} is missing an Orthogroup column")
        species_columns = [column for column in reader.fieldnames if column != "Orthogroup"]
        for row in reader:
            orthogroup_id = row["Orthogroup"]
            orthogroups[orthogroup_id] = {}
            for species_column in species_columns:
                proteins = [value.strip() for value in row.get(species_column, "").split(",") if value.strip()]
                orthogroups[orthogroup_id][species_column] = proteins
                for protein_id in proteins:
                    protein_to_group[protein_id] = (orthogroup_id, species_column)
    return orthogroups, protein_to_group


def result_roots_from_gene_counts(gene_count_long_path: Path) -> list[Path]:
    rows = read_tsv(gene_count_long_path, ORTHOGROUP_GENE_COUNT_LONG_COLUMNS)
    roots = sorted({Path(row["orthofinder_results_dir"]) for row in rows if row.get("orthofinder_results_dir")})
    if not roots:
        raise ValueError(f"{gene_count_long_path} does not record an orthofinder_results_dir")
    return roots


def find_orthogroups_tsv(gene_count_long_path: Path) -> Path:
    roots = result_roots_from_gene_counts(gene_count_long_path)
    candidates: list[Path] = []
    for root in roots:
        direct = root / "Orthogroups" / "Orthogroups.tsv"
        if direct.exists():
            candidates.append(direct)
        candidates.extend(path for path in root.glob("**/Orthogroups.tsv") if path.is_file())
    if not candidates:
        searched = ", ".join(str(root) for root in roots)
        raise FileNotFoundError(f"No Orthogroups.tsv found under parsed OrthoFinder result roots: {searched}")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def build_count_lookup(gene_count_long_path: Path) -> dict[str, dict[str, int]]:
    rows = read_tsv(gene_count_long_path, ORTHOGROUP_GENE_COUNT_LONG_COLUMNS)
    counts: dict[str, dict[str, int]] = {}
    for row in rows:
        orthogroup_id = row["orthogroup_id"]
        counts.setdefault(orthogroup_id, {})[row["orthofinder_species_column"]] = int(row["copy_count"])
    return counts


def manifest_resource_lookup(manifest_path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(manifest_path, ORTHOFINDER_INPUT_MANIFEST_COLUMNS)
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        if row["orthofinder_ready"] != "True":
            continue
        path = Path(row["orthofinder_input_path"])
        for key in {path.name, path.stem, path.name.replace(".faa", ""), path.stem.replace(".faa", "")}:
            lookup[key] = row
    return lookup


def coordinate_lookup(coordinates_path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(coordinates_path, REFERENCE_GENE_COORDINATE_COLUMNS)
    return {
        row["gene_id"]: row
        for row in rows
        if row["resource_id"] == TARGET_RESOURCE_ID and row["parse_status"] == "PARSED_FROM_GFF"
    }


def lookup_protein(
    protein_id: str,
    protein_to_group: dict[str, tuple[str, str]],
) -> tuple[str, str, str]:
    for variant in orthofinder_identifier_variants(protein_id):
        if variant in protein_to_group:
            orthogroup_id, species_column = protein_to_group[variant]
            return orthogroup_id, species_column, variant
    return "NOT_ASSESSED", "NOT_ASSESSED", protein_id


def format_species_counts(
    orthogroup_ids: list[str],
    counts: dict[str, dict[str, int]],
    manifest_lookup: dict[str, dict[str, str]],
) -> str:
    if not orthogroup_ids:
        return "NOT_ASSESSED"
    formatted_groups: list[str] = []
    for orthogroup_id in orthogroup_ids:
        species_counts = counts.get(orthogroup_id, {})
        parts = []
        for species_column, count in sorted(species_counts.items()):
            manifest_row = manifest_lookup.get(species_column, {})
            species_id = manifest_row.get("species_id", species_column)
            parts.append(f"{species_id}={count}")
        formatted_groups.append(f"{orthogroup_id}:{','.join(parts)}")
    return "|".join(formatted_groups)


def supporting_resources(
    orthogroup_ids: list[str],
    counts: dict[str, dict[str, int]],
    manifest_lookup: dict[str, dict[str, str]],
) -> str:
    resources: list[str] = []
    for orthogroup_id in orthogroup_ids:
        for species_column, count in counts.get(orthogroup_id, {}).items():
            if count <= 0:
                continue
            resource_id = manifest_lookup.get(species_column, {}).get("resource_id")
            if resource_id:
                resources.append(resource_id)
    return join_values(resources)


def coordinate_summary(gene_ids: list[str], coordinates: dict[str, dict[str, str]]) -> str:
    parts: list[str] = []
    for gene_id in gene_ids:
        row = coordinates.get(gene_id)
        if not row:
            parts.append(f"{gene_id}:NO_COORDINATES")
            continue
        parts.append(f"{gene_id}:{row['seqid']}:{row['start']}-{row['end']}:{row['strand']}")
    return ";".join(parts) if parts else "NOT_ASSESSED"


def _truth(value: bool) -> str:
    return "yes" if value else "no"


def _artifact_risk(copy_count: int, isoform_risk: bool, coordinate_support: bool, orthology_support: bool) -> str:
    if copy_count == 0:
        return "not_assessable"
    if copy_count > 1 and (isoform_risk or not coordinate_support):
        return "high"
    if not orthology_support or isoform_risk:
        return "moderate"
    return "low"


def _required_validation(mapping_status: str, copy_count: int) -> str:
    if mapping_status == "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH":
        return (
            "ANNOTATION_UNCERTAINTY: do not infer absence; validate with reciprocal similarity searches, "
            "protein-to-genome alignment, and manual annotation review."
        )
    if copy_count > 1:
        return (
            "REQUIRES_VALIDATION: confirm protein-domain integrity, remove isoform inflation, inspect separable "
            "genomic loci, and test cross-resource Greenland shark support before any duplication interpretation."
        )
    return (
        "REQUIRES_VALIDATION: confirm orthology, domain integrity, and annotation consistency before biological "
        "interpretation."
    )


def integrate_real_candidates(
    candidate_panel_path: Path,
    gff_path: Path,
    gene_count_long_path: Path,
    gene_coordinates_path: Path,
    manifest_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    candidates = iter_candidates(read_yaml(candidate_panel_path))
    orthogroups_path = find_orthogroups_tsv(gene_count_long_path)
    result_root = str(orthogroups_path.parent.parent)
    _orthogroups, protein_to_group = parse_orthogroups(orthogroups_path)
    counts = build_count_lookup(gene_count_long_path)
    manifest_lookup = manifest_resource_lookup(manifest_path)
    coordinates = coordinate_lookup(gene_coordinates_path)

    proteins = parse_gff_protein_annotations(gff_path)
    proteins_by_label: dict[str, list[ProteinAnnotation]] = {}
    for protein in proteins:
        proteins_by_label.setdefault(normalise_label(protein.gene_symbol), []).append(protein)

    copy_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    target_species_column = "smic__SMIC_TOKYO_GENOME_2025"

    for candidate in candidates:
        gene_symbol = str(candidate["gene_symbol"]).upper()
        aliases = candidate_aliases(candidate)
        matched: list[ProteinAnnotation] = []
        for alias in aliases:
            matched.extend(proteins_by_label.get(alias, []))
        # Deduplicate while preserving stable protein order.
        matched_by_protein = {protein.protein_id: protein for protein in matched}
        matched = list(matched_by_protein.values())

        protein_mappings = [lookup_protein(protein.protein_id, protein_to_group) for protein in matched]
        orthogroup_ids = sorted({mapping[0] for mapping in protein_mappings if mapping[0] != "NOT_ASSESSED"})
        matched_orthofinder_ids = [mapping[2] for mapping in protein_mappings if mapping[0] != "NOT_ASSESSED"]
        gene_ids = sorted({protein.gene_id for protein in matched if protein.gene_id != "NOT_ASSESSED"})
        copy_count = len(gene_ids)
        species_counts = format_species_counts(orthogroup_ids, counts, manifest_lookup)
        target_protein_count = sum(counts.get(orthogroup_id, {}).get(target_species_column, 0) for orthogroup_id in orthogroup_ids)
        comparator_support = any(
            count > 0 and species_column != target_species_column
            for orthogroup_id in orthogroup_ids
            for species_column, count in counts.get(orthogroup_id, {}).items()
        )
        coordinate_support = bool(gene_ids) and all(gene_id in coordinates for gene_id in gene_ids)
        separable_loci = copy_count > 1 and coordinate_support
        isoform_risk = len(matched) > copy_count or len(orthogroup_ids) > copy_count
        fragmentation_risk = any("fragment" in protein.product.lower() or "partial" in protein.product.lower() for protein in matched)

        if matched and orthogroup_ids:
            mapping_status = "ANNOTATION_SYMBOL_MATCH_ORTHOGROUP_MAPPED"
            annotation_match_level = "exact_gene_symbol_or_alias_from_gff"
        elif matched:
            mapping_status = "ANNOTATION_SYMBOL_MATCH_NO_ORTHOGROUP"
            annotation_match_level = "exact_gene_symbol_or_alias_from_gff"
        else:
            mapping_status = "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH"
            annotation_match_level = "NOT_ASSESSED"

        orthology_support = bool(orthogroup_ids and comparator_support)
        artifact_risk = _artifact_risk(copy_count, isoform_risk, coordinate_support, orthology_support)
        resources = supporting_resources(orthogroup_ids, counts, manifest_lookup)
        if resources == "NOT_ASSESSED" and matched:
            resources = TARGET_RESOURCE_ID

        orthogroup_value = ";".join(orthogroup_ids) if orthogroup_ids else "NOT_ASSESSED"
        gene_id_value = ",".join(gene_ids) if gene_ids else "NOT_ASSESSED"
        protein_id_value = ",".join(matched_orthofinder_ids) if matched_orthofinder_ids else "NOT_ASSESSED"

        copy_rows.append(
            {
                "mechanism": str(candidate["mechanism"]),
                "gene_symbol": gene_symbol,
                "species_id": TARGET_SPECIES_ID,
                "resource_id": TARGET_RESOURCE_ID,
                "orthogroup_id": orthogroup_value,
                "copy_count": str(copy_count),
                "gene_ids": gene_id_value,
                "protein_ids": protein_id_value,
                "orthogroup_target_protein_count": str(target_protein_count) if orthogroup_ids else "NOT_ASSESSED",
                "orthogroup_species_counts": species_counts,
                "annotation_match_level": annotation_match_level,
                "mapping_status": mapping_status,
                "orthofinder_results_dir": result_root,
                "resources_supporting": resources,
                "demo_only": "False",
            }
        )

        audit_rows.append(
            {
                "mechanism": str(candidate["mechanism"]),
                "gene_symbol": gene_symbol,
                "orthogroup_id": orthogroup_value,
                "copy_count": str(copy_count),
                "copy_ids": gene_id_value,
                "protein_ids": protein_id_value,
                "orthology_support": _truth(orthology_support),
                "domain_integrity": "NOT_ASSESSED",
                "separable_loci": _truth(separable_loci),
                "coordinate_support": _truth(coordinate_support),
                "cross_resource_support": "no",
                "isoform_risk": _truth(isoform_risk),
                "fragmentation_risk": _truth(fragmentation_risk),
                "expression_support": "NOT_ASSESSED",
                "annotation_match_level": annotation_match_level,
                "orthogroup_species_counts": species_counts,
                "coordinate_summary": coordinate_summary(gene_ids, coordinates),
                "mapping_status": mapping_status,
                "orthofinder_results_dir": result_root,
                "resources_supporting": resources,
                "artifact_risk": artifact_risk,
                "required_validation": _required_validation(mapping_status, copy_count),
                "demo_only": "False",
            }
        )

    return copy_rows, audit_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Map curated candidates onto real OrthoFinder outputs.")
    parser.add_argument("--candidate-panels", type=Path, required=True)
    parser.add_argument("--annotation-gff", type=Path, required=True)
    parser.add_argument("--orthogroup-gene-counts-long", type=Path, required=True)
    parser.add_argument("--gene-coordinates", type=Path, required=True)
    parser.add_argument("--orthofinder-input-manifest", type=Path, required=True)
    parser.add_argument("--copy-number-output", type=Path, required=True)
    parser.add_argument("--duplication-audit-output", type=Path, required=True)
    args = parser.parse_args()

    configure_logging()
    copy_rows, audit_rows = integrate_real_candidates(
        args.candidate_panels,
        args.annotation_gff,
        args.orthogroup_gene_counts_long,
        args.gene_coordinates,
        args.orthofinder_input_manifest,
    )
    write_tsv(args.copy_number_output, copy_rows, COPY_NUMBER_COLUMNS)
    write_tsv(args.duplication_audit_output, audit_rows, DUPLICATION_AUDIT_COLUMNS)


if __name__ == "__main__":
    main()
