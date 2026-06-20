"""Conservative consolidation of Phase 4 validation outputs."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .evidence import validate_guardrail_language
from .schemas import (
    EVIDENCE_COLUMNS,
    PHASE4C_GENE_REVIEW_COLUMNS,
    PHASE4C_LOCUS_REVIEW_COLUMNS,
    PHASE4C_TP53_SUMMARY_COLUMNS,
    PHASE4D_CANDIDATE_INTERPRETATION_COLUMNS,
)
from .utils import read_tsv, write_tsv

LOGGER = logging.getLogger(__name__)

PHASE4D_GENES = {"H1F0", "FTH1B", "RAD51", "TP53"}

GENE_RELEVANCE = {
    "H1F0": "Plausible chromatin-maintenance lead; hypothesis-generating until orthology and locus context are validated.",
    "FTH1B": "Plausible iron-handling lead; hypothesis-generating until ferritin-family loci are manually resolved.",
    "RAD51": "Plausible DNA repair lead; hypothesis-generating until paralog identity and locus context are validated.",
    "TP53": "NOT_ASSESSED",
}

GENE_NEXT_ACTION = {
    "H1F0": "Manually inspect H1 paralog identity, local coordinates, exon structure, and cross-resource support.",
    "FTH1B": "Manually separate ferritin-family paralogs, isoforms, fragments, and local genomic loci before copy-number interpretation.",
    "RAD51": "Manually resolve RAD51-family paralogy, local coordinates, and cross-resource support.",
    "TP53": "Run additional TP53-family queries, inspect raw miniprot GFF/PAF, and compare against independent resources before interpretation.",
}


def parse_int(value: str | None) -> int:
    if value in {None, "", "NOT_ASSESSED"}:
        return 0
    return int(float(str(value)))


def phase4_sources_for_gene(gene_symbol: str) -> str:
    if gene_symbol == "TP53":
        return "SMIC_TOKYO_GENOME_2025;RTYP_REFSEQ_2022;miniprot_0.18;phase4_forward_hit_target_scaffolds"
    return "SMIC_TOKYO_GENOME_2025;Pfam-A_38.1;miniprot_0.18;phase4b_candidate_scaffolds;phase4c_locus_review"


def phase4_supporting_files(gene_symbol: str) -> str:
    if gene_symbol == "TP53":
        files = [
            "results/rescue/tp53_forward_hit_target_regions.tsv",
            "results/rescue/tp53_targeted_forward_miniprot_raw.gff",
            "results/rescue/tp53_targeted_forward_alignment_hits.tsv",
            "results/rescue/tp53_targeted_forward_search_summary.tsv",
        ]
    else:
        files = [
            "results/rescue/phase4c_locus_review.tsv",
            "results/rescue/phase4c_gene_review_summary.tsv",
            "results/rescue/phase4c_rescue_domain_integrity.tsv",
            "results/rescue/phase4b_genome_alignment_hits.tsv",
            "results/rescue/phase4b_miniprot_raw.gff",
        ]
    return ";".join(files)


def candidate_locus_interpretation(row: dict[str, str]) -> str:
    gene_symbol = row["gene_symbol"]
    if gene_symbol == "FTH1B":
        return (
            "Multiple ferritin-like candidate loci have high-coverage miniprot support and Ferritin Pfam support. "
            "This is candidate-locus evidence only; copy-number expansion is not validated because ferritin-family "
            "paralogy, isoforms, and local assembly context require manual review."
        )
    if gene_symbol == "H1F0":
        return (
            "Multiple H1F0-like candidate loci have high-coverage miniprot support and linker-histone Pfam support. "
            "This is candidate-locus evidence only; H1 paralog identity and duplication status require manual review."
        )
    if gene_symbol == "RAD51":
        return (
            "Multiple RAD51-like candidate loci have high-coverage miniprot support and RAD51-family domain support. "
            "This is candidate-locus evidence only; RAD51 paralog identity and duplication status require manual review."
        )
    return "Candidate-locus evidence requires manual review before biological interpretation."


def build_gene_candidate_interpretation(row: dict[str, str]) -> dict[str, str]:
    gene_symbol = row["gene_symbol"]
    reviewed_loci = parse_int(row.get("reviewed_locus_count"))
    high_coverage_loci = parse_int(row.get("high_coverage_locus_count"))
    domain_supported_loci = parse_int(row.get("domain_supported_locus_count"))
    disruption_loci = parse_int(row.get("loci_with_disruption_flags"))

    if reviewed_loci > 0 and high_coverage_loci > 0 and domain_supported_loci > 0:
        evidence_tier = "Tier 2"
        phase4d_status = "CANDIDATE_LOCI_SUPPORTED_REQUIRES_MANUAL_REVIEW"
        evidence_type = "phase4d_candidate_locus_review"
        copy_status = "NOT_VALIDATED_DUPLICATION"
    else:
        evidence_tier = "Artifact/uncertain"
        phase4d_status = "NOT_ASSESSABLE_FROM_CURRENT_PHASE4_TABLES"
        evidence_type = "phase4d_candidate_locus_review_uncertain"
        copy_status = "NOT_ASSESSED"

    artifact_risk = row.get("artifact_risk", "moderate") or "moderate"
    if gene_symbol == "FTH1B":
        artifact_risk = "high"
    elif disruption_loci > 0 and artifact_risk == "low":
        artifact_risk = "moderate"

    return {
        "gene_symbol": gene_symbol,
        "mechanism": row["mechanism"],
        "phase4d_status": phase4d_status,
        "evidence_type": evidence_type,
        "evidence_tier": evidence_tier,
        "phase4_sources": phase4_sources_for_gene(gene_symbol),
        "candidate_loci": row.get("candidate_loci", "NOT_ASSESSED"),
        "reviewed_locus_count": row.get("reviewed_locus_count", "0"),
        "high_coverage_locus_count": row.get("high_coverage_locus_count", "0"),
        "domain_supported_locus_count": row.get("domain_supported_locus_count", "0"),
        "loci_with_disruption_flags": row.get("loci_with_disruption_flags", "0"),
        "best_query_coverage": row.get("best_query_coverage", "NOT_ASSESSED"),
        "artifact_risk": artifact_risk,
        "conservative_interpretation": candidate_locus_interpretation(row),
        "relevance_to_aging_longevity": GENE_RELEVANCE.get(gene_symbol, "NOT_ASSESSED"),
        "copy_number_interpretation_status": copy_status,
        "required_validation": row.get("required_validation", "REQUIRES_MANUAL_LOCUS_REVIEW_AND_CROSS_RESOURCE_SUPPORT"),
        "claim_language_guardrail": (
            "Do not claim validated duplication, functional advantage, pathway activity, adaptation, or causation from Phase 4d."
        ),
        "supporting_files": phase4_supporting_files(gene_symbol),
        "next_action": GENE_NEXT_ACTION.get(gene_symbol, "Manual review before interpretation."),
        "notes": row.get("notes", ""),
    }


def build_tp53_interpretation(row: dict[str, str]) -> dict[str, str]:
    status = row.get("genome_search_status", "NOT_ASSESSED")
    candidate_loci = row.get("candidate_loci", "NOT_ASSESSED")
    high_coverage = parse_int(row.get("high_coverage_hit_count"))
    partial = parse_int(row.get("partial_hit_count"))

    if status == "POSSIBLE_DISRUPTED_TP53_ALIGNMENT":
        interpretation = (
            "A targeted scaffold search found one high-coverage p53-family alignment with miniprot disruption flags. "
            "This is an uncertain candidate locus requiring manual exon, domain, and orthology review; it does not "
            "establish gene state, function, or a longevity mechanism."
        )
        phase4d_status = "TP53_CANDIDATE_ALIGNMENT_UNCERTAIN_REQUIRES_REVIEW"
    elif high_coverage > 0 or partial > 0:
        interpretation = (
            "A targeted scaffold search found p53-family alignment evidence that remains unresolved under the current filters. "
            "This requires manual review before any gene-level interpretation."
        )
        phase4d_status = "TP53_ALIGNMENT_EVIDENCE_UNCERTAIN"
    else:
        interpretation = (
            "The current targeted scaffold search did not produce interpretable TP53 evidence under the configured filters. "
            "This remains annotation uncertainty and requires broader query/resource validation."
        )
        phase4d_status = "TP53_NOT_ASSESSABLE_UNDER_CURRENT_SEARCH"

    return {
        "gene_symbol": "TP53",
        "mechanism": "p53 pathway",
        "phase4d_status": phase4d_status,
        "evidence_type": "phase4d_targeted_tp53_review",
        "evidence_tier": "Artifact/uncertain",
        "phase4_sources": phase4_sources_for_gene("TP53"),
        "candidate_loci": candidate_loci,
        "reviewed_locus_count": row.get("alignment_hit_count", "0"),
        "high_coverage_locus_count": row.get("high_coverage_hit_count", "0"),
        "domain_supported_locus_count": "NOT_ASSESSED",
        "loci_with_disruption_flags": "1" if status == "POSSIBLE_DISRUPTED_TP53_ALIGNMENT" else "NOT_ASSESSED",
        "best_query_coverage": row.get("best_query_coverage", "NOT_ASSESSED"),
        "artifact_risk": "high",
        "conservative_interpretation": interpretation,
        "relevance_to_aging_longevity": GENE_RELEVANCE["TP53"],
        "copy_number_interpretation_status": "NOT_ASSESSED",
        "required_validation": row.get(
            "required_validation",
            "Use additional TP53-family queries, inspect raw alignments, and compare resources before interpretation.",
        ),
        "claim_language_guardrail": (
            "Do not report TP53 loss, functional state, adaptation, or functional relevance from this targeted search alone."
        ),
        "supporting_files": phase4_supporting_files("TP53"),
        "next_action": GENE_NEXT_ACTION["TP53"],
        "notes": "Phase 4d retains TP53 as a candidate p53-family alignment requiring manual review.",
    }


def build_phase4d_interpretations(
    gene_summary_path: Path,
    locus_review_path: Path,
    tp53_summary_path: Path,
) -> list[dict[str, str]]:
    gene_rows = read_tsv(gene_summary_path, PHASE4C_GENE_REVIEW_COLUMNS)
    locus_rows = read_tsv(locus_review_path, PHASE4C_LOCUS_REVIEW_COLUMNS)
    tp53_rows = read_tsv(tp53_summary_path, PHASE4C_TP53_SUMMARY_COLUMNS)

    reviewed_genes = {row["gene_symbol"] for row in locus_rows}
    interpretations: list[dict[str, str]] = []
    for row in sorted(gene_rows, key=lambda item: item["gene_symbol"]):
        if row["gene_symbol"] not in {"H1F0", "FTH1B", "RAD51"}:
            continue
        if row["gene_symbol"] not in reviewed_genes:
            raise ValueError(f"{row['gene_symbol']} is in the Phase 4c gene summary but not the locus review table")
        interpretations.append(build_gene_candidate_interpretation(row))

    if not tp53_rows:
        raise ValueError(f"{tp53_summary_path} has no TP53 summary rows")
    tp53_row = next((row for row in tp53_rows if row.get("gene_symbol") == "TP53"), tp53_rows[0])
    interpretations.append(build_tp53_interpretation(tp53_row))
    return interpretations


def integrated_row_from_interpretation(row: dict[str, str]) -> dict[str, str]:
    return {
        "mechanism": row["mechanism"],
        "gene_or_pathway": row["gene_symbol"],
        "evidence_type": row["evidence_type"],
        "evidence_tier": row["evidence_tier"],
        "resources_supporting": row["phase4_sources"],
        "artifact_risk": row["artifact_risk"],
        "biological_interpretation": row["conservative_interpretation"],
        "relevance_to_aging_longevity": row["relevance_to_aging_longevity"],
        "translational_relevance": "NOT_ASSESSED",
        "required_validation": row["required_validation"],
        "claim_language_guardrail": row["claim_language_guardrail"],
    }


def consolidate_integrated_evidence(
    base_evidence_path: Path,
    interpretations: list[dict[str, str]],
) -> list[dict[str, str]]:
    base_rows = read_tsv(base_evidence_path, EVIDENCE_COLUMNS)
    retained_rows = [row for row in base_rows if row.get("gene_or_pathway") not in PHASE4D_GENES]
    phase4d_rows = [integrated_row_from_interpretation(row) for row in interpretations]
    consolidated = retained_rows + sorted(phase4d_rows, key=lambda row: row["gene_or_pathway"])
    validate_guardrail_language(consolidated)
    return consolidated


def write_consolidation_outputs(
    base_evidence_path: Path,
    gene_summary_path: Path,
    locus_review_path: Path,
    tp53_summary_path: Path,
    interpretation_output: Path,
    integrated_output: Path,
) -> None:
    interpretations = build_phase4d_interpretations(gene_summary_path, locus_review_path, tp53_summary_path)
    write_tsv(interpretation_output, interpretations, PHASE4D_CANDIDATE_INTERPRETATION_COLUMNS)
    consolidated = consolidate_integrated_evidence(base_evidence_path, interpretations)
    write_tsv(integrated_output, consolidated, EVIDENCE_COLUMNS)
    LOGGER.info("Wrote Phase 4d interpretation for %d genes", len(interpretations))


def main() -> None:
    parser = argparse.ArgumentParser(description="Consolidate Phase 4 validation outputs conservatively.")
    parser.add_argument("--base-evidence", type=Path, required=True)
    parser.add_argument("--phase4c-gene-summary", type=Path, required=True)
    parser.add_argument("--phase4c-locus-review", type=Path, required=True)
    parser.add_argument("--tp53-summary", type=Path, required=True)
    parser.add_argument("--interpretation-output", type=Path, required=True)
    parser.add_argument("--integrated-output", type=Path, required=True)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    write_consolidation_outputs(
        args.base_evidence,
        args.phase4c_gene_summary,
        args.phase4c_locus_review,
        args.tp53_summary,
        args.interpretation_output,
        args.integrated_output,
    )


if __name__ == "__main__":
    main()
