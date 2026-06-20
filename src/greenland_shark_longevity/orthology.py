"""Candidate and toy orthogroup integration."""

from __future__ import annotations

import argparse
import logging
from collections import defaultdict
from pathlib import Path

from .candidate_panels import iter_candidates
from .schemas import COPY_NUMBER_COLUMNS, DUPLICATION_AUDIT_COLUMNS
from .utils import as_bool, join_values, read_tsv, read_yaml, write_tsv


def _split_genes(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def _truth(value: bool) -> str:
    return "yes" if value else "no"


def integrate_candidates(
    candidate_panel_path: Path,
    orthogroups_path: Path,
    gene_coordinates_path: Path,
    domains_path: Path,
    expression_path: Path,
    species_id: str = "smic",
    demo_only: bool = True,
) -> tuple[list[dict], list[dict]]:
    candidates = iter_candidates(read_yaml(candidate_panel_path))
    orthogroups = read_tsv(orthogroups_path, ["orthogroup_id", "species_id", "gene_ids"])
    coordinates = read_tsv(
        gene_coordinates_path,
        ["gene_id", "gene_symbol", "species_id", "resource_id", "scaffold", "start", "end", "locus_id", "isoform_group"],
    )
    domains = read_tsv(domains_path, ["gene_id", "expected_domain_set", "domain_integrity", "domain_notes"])
    expression = read_tsv(expression_path, ["gene_symbol", "species_id", "tissue", "resource_id", "evidence_status", "notes"])

    gene_to_orthogroup: dict[str, str] = {}
    orthogroup_species: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in orthogroups:
        genes = _split_genes(row["gene_ids"])
        orthogroup_species[(row["orthogroup_id"], row["species_id"])].extend(genes)
        for gene_id in genes:
            gene_to_orthogroup[gene_id] = row["orthogroup_id"]

    coords_by_symbol: dict[str, list[dict]] = defaultdict(list)
    for row in coordinates:
        if row["species_id"] == species_id:
            coords_by_symbol[row["gene_symbol"].upper()].append(row)

    domains_by_gene = {row["gene_id"]: row for row in domains}
    expression_by_symbol: dict[str, list[dict]] = defaultdict(list)
    for row in expression:
        if row["species_id"] == species_id:
            expression_by_symbol[row["gene_symbol"].upper()].append(row)

    copy_rows: list[dict] = []
    audit_rows: list[dict] = []
    for candidate in candidates:
        gene_symbol = str(candidate["gene_symbol"]).upper()
        genes = coords_by_symbol.get(gene_symbol, [])
        gene_ids = [row["gene_id"] for row in genes]
        orthogroup_ids = sorted({gene_to_orthogroup.get(gene_id, "NOT_ASSESSED") for gene_id in gene_ids})
        orthogroup_id = ";".join(orthogroup_ids) if orthogroup_ids else "NOT_ASSESSED"
        resources = join_values(row["resource_id"] for row in genes)
        loci = {row["locus_id"] for row in genes if row.get("locus_id")}
        isoforms = {row["isoform_group"] for row in genes if row.get("isoform_group")}
        domain_states = [domains_by_gene.get(gene_id, {}).get("domain_integrity", "NOT_ASSESSED") for gene_id in gene_ids]
        expression_states = [row["evidence_status"] for row in expression_by_symbol.get(gene_symbol, [])]

        orthology_support = bool(gene_ids and orthogroup_id != "NOT_ASSESSED")
        complete_domains = bool(domain_states) and all(state == "complete" for state in domain_states)
        partial_domains = any(state == "partial" for state in domain_states)
        separable_loci = len(loci) >= len(gene_ids) and len(gene_ids) > 1
        isoform_risk = len(gene_ids) > 1 and len(isoforms) < len(gene_ids)
        fragmentation_risk = partial_domains or any("frag" in gene_id.lower() for gene_id in gene_ids)
        artifact_risk = "high" if isoform_risk or fragmentation_risk else ("low" if gene_ids else "not_assessable")

        copy_rows.append(
            {
                "mechanism": candidate["mechanism"],
                "gene_symbol": gene_symbol,
                "species_id": species_id,
                "resource_id": resources,
                "orthogroup_id": orthogroup_id,
                "copy_count": len(gene_ids),
                "gene_ids": ",".join(gene_ids) if gene_ids else "NOT_ASSESSED",
                "protein_ids": "NOT_ASSESSED",
                "orthogroup_target_protein_count": "NOT_ASSESSED",
                "orthogroup_species_counts": "NOT_ASSESSED",
                "annotation_match_level": "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE",
                "mapping_status": "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE",
                "orthofinder_results_dir": "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE",
                "resources_supporting": resources,
                "demo_only": str(demo_only),
            }
        )

        audit_rows.append(
            {
                "mechanism": candidate["mechanism"],
                "gene_symbol": gene_symbol,
                "orthogroup_id": orthogroup_id,
                "copy_count": len(gene_ids),
                "copy_ids": ",".join(gene_ids) if gene_ids else "NOT_ASSESSED",
                "protein_ids": "NOT_ASSESSED",
                "orthology_support": _truth(orthology_support),
                "domain_integrity": "complete" if complete_domains else ("partial" if partial_domains else "NOT_ASSESSED"),
                "separable_loci": _truth(separable_loci),
                "coordinate_support": _truth(bool(gene_ids)),
                "cross_resource_support": "no",
                "isoform_risk": _truth(isoform_risk),
                "fragmentation_risk": _truth(fragmentation_risk),
                "expression_support": ";".join(sorted(set(expression_states))) if expression_states else "NOT_ASSESSED",
                "annotation_match_level": "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE",
                "orthogroup_species_counts": "NOT_ASSESSED",
                "coordinate_summary": "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE",
                "mapping_status": "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE",
                "orthofinder_results_dir": "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE",
                "resources_supporting": resources,
                "artifact_risk": artifact_risk,
                "required_validation": "Replace demo data with verified public resources; confirm orthology, domains, loci, and isoforms.",
                "demo_only": str(demo_only),
            }
        )
    return copy_rows, audit_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Integrate candidate panels with orthogroup and annotation evidence.")
    parser.add_argument("--candidate-panels", type=Path, required=True)
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--gene-coordinates", type=Path, required=True)
    parser.add_argument("--domains", type=Path, required=True)
    parser.add_argument("--expression-support", type=Path, required=True)
    parser.add_argument("--copy-number-output", type=Path, required=True)
    parser.add_argument("--duplication-audit-output", type=Path, required=True)
    parser.add_argument("--species-id", default="smic")
    parser.add_argument("--demo-only", default="true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    copy_rows, audit_rows = integrate_candidates(
        args.candidate_panels,
        args.orthogroups,
        args.gene_coordinates,
        args.domains,
        args.expression_support,
        species_id=args.species_id,
        demo_only=as_bool(args.demo_only),
    )
    write_tsv(args.copy_number_output, copy_rows, COPY_NUMBER_COLUMNS)
    write_tsv(args.duplication_audit_output, audit_rows, DUPLICATION_AUDIT_COLUMNS)


if __name__ == "__main__":
    main()
