"""Phase 4e manual-review hardening for high-priority candidate loci."""

from __future__ import annotations

import argparse
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from .evidence import validate_guardrail_language
from .schemas import (
    EVIDENCE_COLUMNS,
    PHASE4B_ALIGNMENT_HIT_COLUMNS,
    PHASE4C_LOCUS_REVIEW_COLUMNS,
    PHASE4C_TP53_TARGET_REGION_COLUMNS,
    PHASE4E_GENE_HARDENING_COLUMNS,
    PHASE4E_LOCUS_HARDENING_COLUMNS,
)
from .utils import read_tsv, write_tsv

LOGGER = logging.getLogger(__name__)
NOT_ASSESSED = "NOT_ASSESSED"


@dataclass(frozen=True)
class GeneFeature:
    seqid: str
    start: int
    end: int
    strand: str
    gene_id: str
    locus_tag: str
    product: str


def parse_int(value: str | None, default: int = 0) -> int:
    if value in {None, "", NOT_ASSESSED}:
        return default
    return int(float(str(value)))


def parse_float(value: str | None, default: float = 0.0) -> float:
    if value in {None, "", NOT_ASSESSED}:
        return default
    return float(str(value))


def parse_attributes(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in value.split(";"):
        if not part:
            continue
        if "=" in part:
            key, attr_value = part.split("=", 1)
            attrs[key] = unquote(attr_value)
    return attrs


def accession_to_scaffold_alias(seqid: str) -> str | None:
    match = re.fullmatch(r"JBLTJD01(\d+)\.\d+", seqid)
    if not match:
        return None
    return f"scaffold_{int(match.group(1))}"


def split_aliases(value: str) -> list[str]:
    if value in {"", NOT_ASSESSED}:
        return []
    return [part.strip() for part in re.split(r"[,;]", value) if part.strip()]


def annotation_seqids_for_locus(row: dict[str, str]) -> list[str]:
    aliases = split_aliases(row.get("requested_scaffold_aliases", ""))
    converted = accession_to_scaffold_alias(row["target_seqid"])
    if converted:
        aliases.append(converted)
    aliases.append(row["target_seqid"])
    seen: set[str] = set()
    ordered: list[str] = []
    for alias in aliases:
        if alias not in seen:
            ordered.append(alias)
            seen.add(alias)
    return ordered


def locus_tags_from_proteins(value: str) -> set[str]:
    return set(re.findall(r"gs_\d+", value or ""))


def interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    left = max(start_a, start_b)
    right = min(end_a, end_b)
    return max(0, right - left + 1)


def interval_distance(row_a: dict[str, str], row_b: dict[str, str]) -> int | None:
    if row_a["annotation_seqid"] != row_b["annotation_seqid"]:
        return None
    a_start = parse_int(row_a["locus_start"])
    a_end = parse_int(row_a["locus_end"])
    b_start = parse_int(row_b["locus_start"])
    b_end = parse_int(row_b["locus_end"])
    if interval_overlap(a_start, a_end, b_start, b_end):
        return 0
    if a_end < b_start:
        return b_start - a_end - 1
    return a_start - b_end - 1


def load_gene_features(gff_path: Path, seqids: set[str]) -> list[GeneFeature]:
    features: list[GeneFeature] = []
    with gff_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene" or fields[0] not in seqids:
                continue
            attrs = parse_attributes(fields[8])
            features.append(
                GeneFeature(
                    seqid=fields[0],
                    start=parse_int(fields[3]),
                    end=parse_int(fields[4]),
                    strand=fields[6],
                    gene_id=attrs.get("ID", NOT_ASSESSED),
                    locus_tag=attrs.get("locus_tag", attrs.get("Name", NOT_ASSESSED)),
                    product=attrs.get("description", attrs.get("product", NOT_ASSESSED)),
                )
            )
    return features


def overlapping_features(
    features: list[GeneFeature],
    seqid_candidates: list[str],
    start: int,
    end: int,
) -> tuple[str, list[GeneFeature], float]:
    locus_span = max(1, end - start + 1)
    for seqid in seqid_candidates:
        overlaps = [feature for feature in features if feature.seqid == seqid and interval_overlap(start, end, feature.start, feature.end)]
        if overlaps:
            max_overlap = max(interval_overlap(start, end, feature.start, feature.end) for feature in overlaps) / locus_span
            return seqid, overlaps, max_overlap
    return (seqid_candidates[0] if seqid_candidates else NOT_ASSESSED, [], 0.0)


def gene_family_risk(gene_symbol: str) -> str:
    risks = {
        "FTH1B": "high_ferritin_family_cluster",
        "H1F0": "moderate_histone_paralogy",
        "RAD51": "moderate_rad51_paralogy",
        "TP53": "high_p53_family_paralogy",
    }
    return risks.get(gene_symbol, "not_assessed")


def product_consistency(gene_symbol: str, products: list[str]) -> str:
    if not products:
        return "NO_ANNOTATION_OVERLAP"
    text = " ".join(products).lower()
    patterns = {
        "FTH1B": ["ferritin"],
        "H1F0": ["histone h1", "histone h5", "linker histone"],
        "RAD51": ["rad51"],
        "TP53": ["p53", "tumor antigen"],
    }
    if any(pattern in text for pattern in patterns.get(gene_symbol, [])):
        return "PRODUCT_CONSISTENT"
    return "PRODUCT_NOT_CONSISTENT"


def base_artifact_flags(
    gene_symbol: str,
    disruption: bool,
    focal_overlap: bool,
    product_status: str,
    overlaps: list[GeneFeature],
    nearest_distance: str,
) -> list[str]:
    flags: list[str] = []
    if gene_symbol == "FTH1B":
        flags.append("FERRITIN_FAMILY_PARALOGY")
        flags.append("CLUSTERED_CANDIDATE_LOCI")
    elif gene_symbol == "H1F0":
        flags.append("HISTONE_PARALOGY")
    elif gene_symbol == "RAD51":
        flags.append("RAD51_FAMILY_PARALOGY")
    elif gene_symbol == "TP53":
        flags.append("P53_FAMILY_PARALOGY")
    if disruption:
        flags.append("MINIPROT_DISRUPTION_TAG")
    if not overlaps:
        flags.append("NO_OVERLAPPING_EXISTING_GENE_MODEL")
    elif not focal_overlap and product_status == "PRODUCT_CONSISTENT":
        flags.append("PRODUCT_CONSISTENT_NONFOCAL_ANNOTATION")
    elif not focal_overlap:
        flags.append("OVERLAPS_NONFOCAL_ANNOTATION")
    if nearest_distance not in {"", NOT_ASSESSED} and parse_int(nearest_distance) < 10000:
        flags.append("NEARBY_SAME_GENE_CANDIDATE_LOCUS")
    return flags or ["NO_SPECIFIC_FLAG_BEYOND_MANUAL_REVIEW"]


def locus_status(
    gene_symbol: str,
    domain_status: str,
    disruption: bool,
    focal_overlap: bool,
    product_status: str,
    overlaps: list[GeneFeature],
) -> str:
    if gene_symbol == "TP53":
        return "UNRESOLVED_P53_FAMILY_ALIGNMENT_REQUIRES_MANUAL_REVIEW"
    if disruption:
        return "POSSIBLE_ALIGNMENT_DISRUPTION_REQUIRES_MANUAL_REVIEW"
    if domain_status == "DOMAIN_SUPPORTED" and focal_overlap:
        return "ANNOTATION_AND_DOMAIN_SUPPORTED_CANDIDATE_LOCUS"
    if domain_status == "DOMAIN_SUPPORTED" and product_status == "PRODUCT_CONSISTENT":
        return "DOMAIN_SUPPORTED_PRODUCT_CONSISTENT_CANDIDATE_LOCUS"
    if domain_status == "DOMAIN_SUPPORTED" and overlaps:
        return "DOMAIN_SUPPORTED_NONFOCAL_ANNOTATION_OVERLAP"
    if domain_status == "DOMAIN_SUPPORTED":
        return "DOMAIN_SUPPORTED_ALIGNMENT_WITHOUT_EXISTING_GENE_OVERLAP"
    return "LOCUS_REQUIRES_MANUAL_REVIEW"


def locus_interpretation(gene_symbol: str, focal_overlap: bool, disruption: bool) -> str:
    if gene_symbol == "TP53":
        return (
            "Targeted p53-family miniprot alignment has domain support but remains unresolved because identity, exon structure, "
            "and parsed disruption flags require manual review."
        )
    if gene_symbol == "FTH1B":
        return (
            "Ferritin-like candidate locus is supported by miniprot and Pfam; it remains a candidate because clustered "
            "ferritin-family paralogy and local annotation context must be manually resolved."
        )
    if disruption:
        return "Candidate locus has miniprot support but parsed disruption flags prevent stronger interpretation."
    if focal_overlap:
        return "Candidate locus overlaps the same annotation model that supplied the rescue protein and has domain support; duplication status remains unvalidated."
    return "Candidate locus has validation support but annotation overlap is incomplete or nonfocal; manual review remains required."


def review_priority(gene_symbol: str, disruption: bool, focal_overlap: bool) -> str:
    if gene_symbol in {"TP53", "FTH1B"} or disruption or not focal_overlap:
        return "HIGH"
    return "MEDIUM"


def make_tp53_locus_rows(alignment_path: Path, target_regions_path: Path, domain_by_query: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    target_alias_by_accession: dict[str, str] = {}
    if target_regions_path.exists():
        for row in read_tsv(target_regions_path, PHASE4C_TP53_TARGET_REGION_COLUMNS):
            header = row.get("target_fasta_header", "")
            if header and header != NOT_ASSESSED:
                target_alias_by_accession[header.split()[0]] = row.get("requested_seqid", NOT_ASSESSED)

    rows: list[dict[str, str]] = []
    for idx, hit in enumerate(read_tsv(alignment_path, PHASE4B_ALIGNMENT_HIT_COLUMNS), start=1):
        query_id = hit["query_id"]
        domain = domain_by_query.get(query_id, {})
        rows.append(
            {
                "gene_symbol": "TP53",
                "mechanism": "p53 pathway",
                "locus_cluster_id": f"TP53_L{idx:03d}",
                "target_seqid": hit["target_seqid"],
                "requested_scaffold_aliases": target_alias_by_accession.get(hit["target_seqid"], NOT_ASSESSED),
                "cluster_start": str(min(parse_int(hit["target_start"]), parse_int(hit["target_end"]))),
                "cluster_end": str(max(parse_int(hit["target_start"]), parse_int(hit["target_end"]))),
                "strand_values": hit["strand"],
                "supporting_query_count": "1",
                "supporting_query_ids": query_id,
                "supporting_original_protein_ids": hit["original_protein_id"],
                "best_query_coverage": hit["query_coverage"],
                "best_identity": hit["identity"],
                "cds_feature_count_range": hit["cds_feature_count"],
                "paf_frameshift_count": "1" if hit["frameshift_or_stop_flag"] == "True" else "0",
                "paf_stop_count": "0",
                "miniprot_disruption_status": "POTENTIAL_DISRUPTION_TAG" if hit["frameshift_or_stop_flag"] == "True" else "NO_PARSED_FS_OR_ST_TAG",
                "domain_validation_status": domain.get("domain_validation_status", NOT_ASSESSED),
                "best_domain_hits": f"{query_id}:{domain.get('domain_validation_status', NOT_ASSESSED)}:{domain.get('best_pfam_accession', NOT_ASSESSED)}:{domain.get('best_pfam_name', NOT_ASSESSED)}",
            }
        )
    return rows


def supporting_files_for_gene(gene_symbol: str) -> str:
    files = [
        "results/rescue/phase4e_locus_manual_review.tsv",
        "results/rescue/phase4e_gene_hardened_summary.tsv",
        "results/rescue/phase4c_locus_review.tsv",
        "results/rescue/phase4c_rescue_domain_integrity.tsv",
    ]
    if gene_symbol == "TP53":
        files.extend(
            [
                "results/rescue/tp53_targeted_forward_alignment_hits.tsv",
                "results/rescue/tp53_targeted_forward_miniprot_raw.gff",
            ]
        )
    else:
        files.append("results/rescue/phase4b_miniprot_raw.gff")
    return ";".join(files)


def build_phase4e_locus_rows(
    phase4c_locus_path: Path,
    domain_integrity_path: Path,
    tp53_alignment_path: Path,
    tp53_target_regions_path: Path,
    annotation_gff_path: Path,
) -> list[dict[str, str]]:
    domain_by_query = {}
    if domain_integrity_path.exists():
        for row in read_tsv(domain_integrity_path):
            domain_by_query[row["representative_fasta_id"]] = row

    base_rows = read_tsv(phase4c_locus_path, PHASE4C_LOCUS_REVIEW_COLUMNS)
    work_rows = base_rows + make_tp53_locus_rows(tp53_alignment_path, tp53_target_regions_path, domain_by_query)

    needed_seqids: set[str] = set()
    for row in work_rows:
        needed_seqids.update(annotation_seqids_for_locus(row))
    gene_features = load_gene_features(annotation_gff_path, needed_seqids)

    preliminary_rows: list[dict[str, str]] = []
    for row in work_rows:
        start = parse_int(row["cluster_start"])
        end = parse_int(row["cluster_end"])
        seqid_candidates = annotation_seqids_for_locus(row)
        annotation_seqid, overlaps, overlap_fraction = overlapping_features(gene_features, seqid_candidates, start, end)
        expected_tags = locus_tags_from_proteins(row.get("supporting_original_protein_ids", ""))
        overlapping_tags = {feature.locus_tag for feature in overlaps}
        overlapping_products = sorted({feature.product for feature in overlaps})
        product_status = product_consistency(row["gene_symbol"], overlapping_products)
        focal_overlap = bool(expected_tags & overlapping_tags)
        disruption = row.get("miniprot_disruption_status") == "POTENTIAL_DISRUPTION_TAG" or parse_int(row.get("paf_frameshift_count")) > 0 or parse_int(row.get("paf_stop_count")) > 0
        locus_span = max(1, end - start + 1)
        overlap_status = (
            "OVERLAPS_FOCAL_ANNOTATION_MODEL"
            if focal_overlap
            else ("OVERLAPS_NONFOCAL_ANNOTATION" if overlaps else "NO_OVERLAPPING_ANNOTATED_GENE")
        )
        preliminary_rows.append(
            {
                "gene_symbol": row["gene_symbol"],
                "mechanism": row["mechanism"],
                "locus_cluster_id": row["locus_cluster_id"],
                "target_seqid": row["target_seqid"],
                "annotation_seqid": annotation_seqid,
                "locus_start": str(start),
                "locus_end": str(end),
                "locus_span_bp": str(locus_span),
                "strand_values": row["strand_values"],
                "same_gene_locus_order": "0",
                "nearest_same_gene_locus_distance_bp": NOT_ASSESSED,
                "supporting_query_count": row["supporting_query_count"],
                "supporting_query_ids": row["supporting_query_ids"],
                "supporting_original_protein_ids": row["supporting_original_protein_ids"],
                "best_query_coverage": row["best_query_coverage"],
                "best_identity": row["best_identity"],
                "cds_feature_count_range": row["cds_feature_count_range"],
                "miniprot_frameshift_count": row.get("paf_frameshift_count", "0"),
                "miniprot_stop_count": row.get("paf_stop_count", "0"),
                "miniprot_disruption_status": "POTENTIAL_DISRUPTION_TAG" if disruption else "NO_PARSED_FS_OR_ST_TAG",
                "domain_validation_status": row["domain_validation_status"],
                "best_domain_hits": row["best_domain_hits"],
                "annotation_overlap_status": overlap_status,
                "overlapping_gene_count": str(len(overlaps)),
                "overlapping_gene_ids": ";".join(sorted({feature.gene_id for feature in overlaps})) if overlaps else NOT_ASSESSED,
                "overlapping_locus_tags": ";".join(sorted({feature.locus_tag for feature in overlaps})) if overlaps else NOT_ASSESSED,
                "overlapping_products": ";".join(overlapping_products) if overlaps else NOT_ASSESSED,
                "focal_annotation_overlap": "yes" if focal_overlap else "no",
                "annotation_product_consistency": product_status,
                "annotation_overlap_fraction": f"{overlap_fraction:.6g}",
                "gene_family_artifact_risk": gene_family_risk(row["gene_symbol"]),
                "phase4e_locus_status": locus_status(row["gene_symbol"], row["domain_validation_status"], disruption, focal_overlap, product_status, overlaps),
                "artifact_flags": "",
                "hardened_interpretation": locus_interpretation(row["gene_symbol"], focal_overlap, disruption),
                "manual_review_priority": review_priority(row["gene_symbol"], disruption, focal_overlap),
                "required_validation": "Manual GFF/PAF and CDS/exon review, paralog identity checks, local context/synteny, and cross-resource validation.",
                "claim_language_guardrail": "Do not claim validated duplication, function, adaptation, activation, inactivation, loss, or causation from Phase 4e alone.",
                "supporting_files": supporting_files_for_gene(row["gene_symbol"]),
            }
        )

    rows_by_gene: dict[str, list[dict[str, str]]] = {}
    for row in preliminary_rows:
        rows_by_gene.setdefault(row["gene_symbol"], []).append(row)
    for gene_rows in rows_by_gene.values():
        gene_rows.sort(key=lambda item: (item["annotation_seqid"], parse_int(item["locus_start"]), parse_int(item["locus_end"])))
        for idx, row in enumerate(gene_rows, start=1):
            distances = [
                distance
                for other in gene_rows
                if other is not row
                for distance in [interval_distance(row, other)]
                if distance is not None
            ]
            nearest = str(min(distances)) if distances else NOT_ASSESSED
            row["same_gene_locus_order"] = str(idx)
            row["nearest_same_gene_locus_distance_bp"] = nearest
            row["artifact_flags"] = ";".join(
                base_artifact_flags(
                    row["gene_symbol"],
                    row["miniprot_disruption_status"] == "POTENTIAL_DISRUPTION_TAG",
                    row["focal_annotation_overlap"] == "yes",
                    row["annotation_product_consistency"],
                    [] if row["overlapping_gene_count"] == "0" else [GeneFeature(row["annotation_seqid"], 0, 0, ".", "", "", "")],
                    nearest,
                )
            )
    return sorted(preliminary_rows, key=lambda item: (item["gene_symbol"], parse_int(item["same_gene_locus_order"])))


def summarize_gene_rows(locus_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows_by_gene: dict[str, list[dict[str, str]]] = {}
    for row in locus_rows:
        rows_by_gene.setdefault(row["gene_symbol"], []).append(row)

    summaries: list[dict[str, str]] = []
    for gene_symbol, rows in sorted(rows_by_gene.items()):
        focal_count = sum(1 for row in rows if row["focal_annotation_overlap"] == "yes")
        product_count = sum(1 for row in rows if row["annotation_product_consistency"] == "PRODUCT_CONSISTENT")
        domain_count = sum(1 for row in rows if row["domain_validation_status"] == "DOMAIN_SUPPORTED")
        disruption_count = sum(1 for row in rows if row["miniprot_disruption_status"] == "POTENTIAL_DISRUPTION_TAG")
        no_disruption = len(rows) - disruption_count
        seqids = sorted({row["annotation_seqid"] for row in rows if row["annotation_seqid"] != NOT_ASSESSED})
        distances = [parse_int(row["nearest_same_gene_locus_distance_bp"]) for row in rows if row["nearest_same_gene_locus_distance_bp"] not in {"", NOT_ASSESSED}]
        min_distance = min(distances) if distances else None
        clustered = "yes" if min_distance is not None and min_distance < 100000 else ("no" if len(rows) > 1 else NOT_ASSESSED)

        if gene_symbol == "TP53":
            status = "UNRESOLVED_P53_FAMILY_CANDIDATE_REQUIRES_ADDITIONAL_PHASE4_VALIDATION"
            tier = "Artifact/uncertain"
            artifact_risk = "high"
            copy_status = "NOT_ASSESSED"
            interpretation = (
                "TP53 remains unresolved: the targeted p53-family alignment has domain support but high artifact risk due to "
                "p53-family paralogy, product-consistent but nonfocal annotation context, and miniprot disruption flags."
            )
            next_action = "Run additional p53-family queries and manually inspect raw GFF/PAF, domain coordinates, and cross-resource support."
        elif gene_symbol == "FTH1B":
            status = "HARDENED_CANDIDATE_CLUSTER_REQUIRES_FERRITIN_FAMILY_RESOLUTION"
            tier = "Tier 2"
            artifact_risk = "high"
            copy_status = "NOT_VALIDATED_DUPLICATION"
            interpretation = (
                "FTH1B-like ferritin candidate loci are supported by miniprot, Pfam, and ferritin-family annotation context, but clustered "
                "paralogous context prevents a validated copy-number interpretation."
            )
            next_action = "Resolve ferritin-family paralogs, inspect local repeats/context, and compare independent Greenland shark resources."
        elif focal_count == len(rows) and domain_count == len(rows) and disruption_count == 0:
            status = "HARDENED_CANDIDATE_LOCI_REQUIRE_CROSS_RESOURCE_VALIDATION"
            tier = "Tier 2"
            artifact_risk = "moderate"
            copy_status = "NOT_VALIDATED_DUPLICATION"
            interpretation = (
                f"{gene_symbol} candidate loci are supported by miniprot, Pfam, and focal annotation overlap, but remain candidate "
                "loci until paralog identity, local context, and cross-resource consistency are validated."
            )
            next_action = "Inspect paralog identity, local context/synteny, and cross-resource support before Phase 8 scoring."
        else:
            status = "CANDIDATE_SIGNAL_REQUIRES_MANUAL_REVIEW"
            tier = "Artifact/uncertain"
            artifact_risk = "high"
            copy_status = "NOT_ASSESSED"
            interpretation = "The current locus evidence is incomplete or artifact-prone and requires manual review before interpretation."
            next_action = "Resolve annotation overlap, domain evidence, and local context before scoring."

        summaries.append(
            {
                "gene_symbol": gene_symbol,
                "mechanism": rows[0]["mechanism"],
                "reviewed_locus_count": str(len(rows)),
                "focal_annotation_overlap_count": str(focal_count),
                "product_consistent_annotation_overlap_count": str(product_count),
                "domain_supported_locus_count": str(domain_count),
                "no_disruption_locus_count": str(no_disruption),
                "disruption_flag_locus_count": str(disruption_count),
                "distinct_annotation_seqids": ";".join(seqids) if seqids else NOT_ASSESSED,
                "min_inter_locus_distance_bp": str(min_distance) if min_distance is not None else NOT_ASSESSED,
                "clustered_same_scaffold": clustered,
                "phase4e_hardened_status": status,
                "evidence_tier_recommendation": tier,
                "artifact_risk": artifact_risk,
                "copy_number_interpretation_status": copy_status,
                "conservative_interpretation": interpretation,
                "next_action": next_action,
                "required_validation": "Manual locus review, paralog resolution, local genomic context, and cross-resource validation before biological claims.",
                "claim_language_guardrail": "Do not claim validated duplication, function, adaptation, activation, inactivation, loss, causation, or translational relevance from Phase 4e alone.",
                "supporting_files": supporting_files_for_gene(gene_symbol),
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
            for row in summaries
        ]
    )
    return summaries


def evidence_rows_from_summary(summary_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    evidence_rows: list[dict[str, str]] = []
    for row in summary_rows:
        evidence_rows.append(
            {
                "mechanism": row["mechanism"],
                "gene_or_pathway": row["gene_symbol"],
                "evidence_type": "phase4e_manual_locus_hardening",
                "evidence_tier": row["evidence_tier_recommendation"],
                "resources_supporting": "SMIC_TOKYO_GENOME_2025;Pfam-A_38.1;miniprot_0.18;figshare_annotation",
                "artifact_risk": row["artifact_risk"],
                "biological_interpretation": row["conservative_interpretation"],
                "relevance_to_aging_longevity": "Hypothesis-generating only; requires cross-resource and functional validation." if row["evidence_tier_recommendation"] == "Tier 2" else NOT_ASSESSED,
                "translational_relevance": NOT_ASSESSED,
                "required_validation": row["required_validation"],
                "claim_language_guardrail": row["claim_language_guardrail"],
            }
        )
    validate_guardrail_language(evidence_rows)
    return evidence_rows


def write_phase4e_outputs(
    phase4c_locus: Path,
    domain_integrity: Path,
    tp53_alignment: Path,
    tp53_target_regions: Path,
    annotation_gff: Path,
    locus_output: Path,
    summary_output: Path,
    evidence_output: Path,
) -> None:
    locus_rows = build_phase4e_locus_rows(phase4c_locus, domain_integrity, tp53_alignment, tp53_target_regions, annotation_gff)
    summary_rows = summarize_gene_rows(locus_rows)
    write_tsv(locus_output, locus_rows, PHASE4E_LOCUS_HARDENING_COLUMNS)
    write_tsv(summary_output, summary_rows, PHASE4E_GENE_HARDENING_COLUMNS)
    write_tsv(evidence_output, evidence_rows_from_summary(summary_rows), EVIDENCE_COLUMNS)
    LOGGER.info("Wrote Phase 4e hardening for %d loci and %d genes", len(locus_rows), len(summary_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 4e manual-review hardening for high-priority loci.")
    parser.add_argument("--phase4c-locus-review", type=Path, required=True)
    parser.add_argument("--domain-integrity", type=Path, required=True)
    parser.add_argument("--tp53-alignment-hits", type=Path, required=True)
    parser.add_argument("--tp53-target-regions", type=Path, required=True)
    parser.add_argument("--annotation-gff", type=Path, required=True)
    parser.add_argument("--locus-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--evidence-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    write_phase4e_outputs(
        args.phase4c_locus_review,
        args.domain_integrity,
        args.tp53_alignment_hits,
        args.tp53_target_regions,
        args.annotation_gff,
        args.locus_output,
        args.summary_output,
        args.evidence_output,
    )


if __name__ == "__main__":
    main()
