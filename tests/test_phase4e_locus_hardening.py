from pathlib import Path

from greenland_shark_longevity.phase4e_locus_hardening import (
    accession_to_scaffold_alias,
    build_phase4e_locus_rows,
    summarize_gene_rows,
)
from greenland_shark_longevity.schemas import (
    DOMAIN_INTEGRITY_COLUMNS,
    PHASE4B_ALIGNMENT_HIT_COLUMNS,
    PHASE4C_LOCUS_REVIEW_COLUMNS,
    PHASE4C_TP53_TARGET_REGION_COLUMNS,
)
from greenland_shark_longevity.utils import write_tsv


def test_accession_to_scaffold_alias_uses_last_assembly_number():
    assert accession_to_scaffold_alias("JBLTJD010000033.1") == "scaffold_33"
    assert accession_to_scaffold_alias("JBLTJD010000247.1") == "scaffold_247"
    assert accession_to_scaffold_alias("scaffold_33") is None


def write_minimal_phase4e_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path]:
    phase4c_locus = tmp_path / "phase4c_locus.tsv"
    domains = tmp_path / "domains.tsv"
    tp53_hits = tmp_path / "tp53_hits.tsv"
    tp53_regions = tmp_path / "tp53_regions.tsv"
    gff = tmp_path / "annotation.gff"

    write_tsv(
        phase4c_locus,
        [
            {
                "gene_symbol": "FTH1B",
                "mechanism": "Ferroptosis, iron handling, and oxidative stress",
                "locus_cluster_id": "FTH1B_L001",
                "target_seqid": "JBLTJD010000033.1",
                "requested_scaffold_aliases": "scaffold_33",
                "cluster_start": "100",
                "cluster_end": "200",
                "strand_values": "+",
                "supporting_query_count": "1",
                "supporting_query_ids": "FTH1B__q1",
                "supporting_original_protein_ids": "gnl|WGS:ZZZZ|gs_018941-P1",
                "best_query_coverage": "1",
                "best_identity": "0.99",
                "hit_count": "1",
                "high_coverage_hit_count": "1",
                "cds_feature_count_range": "4",
                "exon_feature_count_range": "0",
                "paf_frameshift_count": "0",
                "paf_stop_count": "0",
                "miniprot_disruption_status": "NO_PARSED_FS_OR_ST_TAG",
                "domain_validation_status": "DOMAIN_SUPPORTED",
                "domain_supported_query_count": "1",
                "best_domain_hits": "FTH1B__q1:DOMAIN_SUPPORTED:PF00210:Ferritin",
                "overlap_cluster_size": "1",
                "overlap_note": "test",
                "candidate_locus_status": "HIGH_COVERAGE_DOMAIN_SUPPORTED_CANDIDATE_LOCUS",
                "artifact_risk": "moderate",
                "claim_language_guardrail": "test",
                "required_validation": "test",
                "notes": "test",
            }
        ],
        PHASE4C_LOCUS_REVIEW_COLUMNS,
    )
    write_tsv(
        domains,
        [
            {
                "mechanism": "p53 pathway",
                "gene_symbol": "TP53",
                "resource_id": "RTYP_REFSEQ_2022",
                "gene_id": "XP_1",
                "orthogroup_id": "NOT_ASSESSED",
                "representative_protein_id": "XP_1",
                "representative_fasta_id": "TP53__q1",
                "representative_length_aa": "357",
                "domain_validation_status": "DOMAIN_SUPPORTED",
                "best_pfam_accession": "PF00870",
                "best_pfam_name": "P53",
            }
        ],
        DOMAIN_INTEGRITY_COLUMNS,
    )
    write_tsv(
        tp53_hits,
        [
            {
                "gene_symbol": "TP53",
                "query_id": "TP53__q1",
                "query_source_type": "comparator_query_for_unresolved_gene",
                "original_protein_id": "XP_1",
                "source_species_id": "rtyp",
                "source_resource_id": "RTYP_REFSEQ_2022",
                "target_seqid": "JBLTJD010000247.1",
                "target_start": "500",
                "target_end": "900",
                "strand": "-",
                "feature_id": "MP1",
                "feature_type": "mRNA",
                "query_start": "1",
                "query_end": "300",
                "query_length_aa": "357",
                "query_coverage": "0.84",
                "cds_feature_count": "6",
                "exon_feature_count": "0",
                "frameshift_or_stop_flag": "True",
                "identity": "0.58",
                "positive": "0.65",
                "score": "843",
                "rank": "1",
                "alignment_status": "POSSIBLE_DISRUPTION",
                "raw_attributes": "Frameshift=1",
                "notes": "test",
            }
        ],
        PHASE4B_ALIGNMENT_HIT_COLUMNS,
    )
    write_tsv(
        tp53_regions,
        [
            {
                "gene_symbol": "TP53",
                "source_table": "test",
                "source_protein_id": "p53",
                "source_description": "p53-like",
                "source_coordinate": "gene:scaffold_247:500-900:-",
                "requested_seqid": "scaffold_247",
                "target_fasta_header": "JBLTJD010000247.1 test scaffold_247",
                "target_fasta_path": "target.fna",
                "extraction_status": "EXTRACTED",
                "notes": "test",
            }
        ],
        PHASE4C_TP53_TARGET_REGION_COLUMNS,
    )
    gff.write_text(
        "##gff-version 3\n"
        "scaffold_33\tGnomon\tgene\t90\t210\t.\t+\t.\tID=gene-gs_018941;Name=gs_018941;description=ferritin heavy chain b-like;gbkey=Gene;gene_biotype=protein_coding;locus_tag=gs_018941\n"
        "scaffold_247\tGnomon\tgene\t450\t950\t.\t-\t.\tID=gene-gs_999999;Name=gs_999999;description=nonfocal p53-family-like region;gbkey=Gene;gene_biotype=protein_coding;locus_tag=gs_999999\n",
        encoding="utf-8",
    )
    return phase4c_locus, domains, tp53_hits, tp53_regions, gff


def test_phase4e_detects_focal_annotation_overlap_and_tp53_uncertainty(tmp_path: Path):
    phase4c_locus, domains, tp53_hits, tp53_regions, gff = write_minimal_phase4e_inputs(tmp_path)

    locus_rows = build_phase4e_locus_rows(phase4c_locus, domains, tp53_hits, tp53_regions, gff)
    fth1b = next(row for row in locus_rows if row["gene_symbol"] == "FTH1B")
    tp53 = next(row for row in locus_rows if row["gene_symbol"] == "TP53")

    assert fth1b["annotation_overlap_status"] == "OVERLAPS_FOCAL_ANNOTATION_MODEL"
    assert fth1b["annotation_product_consistency"] == "PRODUCT_CONSISTENT"
    assert fth1b["phase4e_locus_status"] == "ANNOTATION_AND_DOMAIN_SUPPORTED_CANDIDATE_LOCUS"
    assert "FERRITIN_FAMILY_PARALOGY" in fth1b["artifact_flags"]
    assert tp53["phase4e_locus_status"] == "UNRESOLVED_P53_FAMILY_ALIGNMENT_REQUIRES_MANUAL_REVIEW"
    assert tp53["focal_annotation_overlap"] == "no"
    assert tp53["annotation_product_consistency"] == "PRODUCT_CONSISTENT"
    assert "MINIPROT_DISRUPTION_TAG" in tp53["artifact_flags"]


def test_phase4e_gene_summary_keeps_copy_number_unvalidated(tmp_path: Path):
    phase4c_locus, domains, tp53_hits, tp53_regions, gff = write_minimal_phase4e_inputs(tmp_path)
    locus_rows = build_phase4e_locus_rows(phase4c_locus, domains, tp53_hits, tp53_regions, gff)

    summary = summarize_gene_rows(locus_rows)
    fth1b = next(row for row in summary if row["gene_symbol"] == "FTH1B")
    tp53 = next(row for row in summary if row["gene_symbol"] == "TP53")

    assert fth1b["evidence_tier_recommendation"] == "Tier 2"
    assert fth1b["product_consistent_annotation_overlap_count"] == "1"
    assert fth1b["copy_number_interpretation_status"] == "NOT_VALIDATED_DUPLICATION"
    assert tp53["evidence_tier_recommendation"] == "Artifact/uncertain"
    assert tp53["copy_number_interpretation_status"] == "NOT_ASSESSED"
