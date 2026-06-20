from pathlib import Path

from greenland_shark_longevity.phase4c_locus_review import (
    build_rescue_domain_targets,
    cluster_hits,
    review_loci,
    tp53_forward_target_scaffolds,
    tp53_summary_row,
)
from greenland_shark_longevity.schemas import (
    DOMAIN_INTEGRITY_COLUMNS,
    PHASE4_RECIPROCAL_HIT_COLUMNS,
    PHASE4B_ALIGNMENT_HIT_COLUMNS,
    PHASE4B_QUERY_COLUMNS,
    PHASE4B_TARGET_REGION_COLUMNS,
)
from greenland_shark_longevity.utils import write_tsv


def test_rescue_domain_targets_use_phase4b_query_ids(tmp_path: Path):
    query_table = tmp_path / "queries.tsv"
    query_fasta = tmp_path / "queries.faa"
    write_tsv(
        query_table,
        [
            {
                "gene_symbol": "TP53",
                "mechanism": "p53 pathway",
                "query_id": "TP53__comparator_query__XP_1",
                "query_source_type": "comparator_query_for_unresolved_gene",
                "original_protein_id": "XP_1",
                "original_description": "cellular tumor antigen p53-like",
                "source_species_id": "rtyp",
                "source_resource_id": "RTYP_REFSEQ_2022",
                "source_fasta": "rtyp.faa",
                "sequence_length_aa": "357",
                "selection_reason": "test",
                "query_status": "READY_FOR_GENOME_ALIGNMENT",
                "notes": "test",
            }
        ],
        PHASE4B_QUERY_COLUMNS,
    )

    rows = build_rescue_domain_targets(query_table, query_fasta)
    assert rows[0]["representative_fasta_id"] == "TP53__comparator_query__XP_1"
    assert rows[0]["domain_check_status"] == "READY_FOR_DOMAIN_SCAN"


def test_cluster_hits_groups_overlapping_intervals():
    hits = [
        {"gene_symbol": "FTH1B", "target_seqid": "scaf1", "target_start": "10", "target_end": "100", "query_id": "q1"},
        {"gene_symbol": "FTH1B", "target_seqid": "scaf1", "target_start": "20", "target_end": "90", "query_id": "q2"},
        {"gene_symbol": "FTH1B", "target_seqid": "scaf1", "target_start": "500", "target_end": "600", "query_id": "q3"},
    ]
    clusters = cluster_hits(hits, min_overlap_fraction=0.50)
    assert len(clusters) == 2
    assert sorted(len(cluster) for cluster in clusters) == [1, 2]


def test_review_loci_uses_domain_support_and_paf_disruption_tags(tmp_path: Path):
    config = tmp_path / "config.yaml"
    hits = tmp_path / "hits.tsv"
    raw_gff = tmp_path / "raw.gff"
    domains = tmp_path / "domains.tsv"
    target_regions = tmp_path / "target_regions.tsv"
    config.write_text(
        "phase4c_locus_review:\n"
        "  review_genes: [FTH1B]\n"
        "  min_overlap_fraction_for_same_locus: 0.50\n",
        encoding="utf-8",
    )
    write_tsv(
        hits,
        [
            {
                "gene_symbol": "FTH1B",
                "query_id": "FTH1B__q1",
                "query_source_type": "focal_rescue_candidate",
                "original_protein_id": "prot1",
                "source_species_id": "smic",
                "source_resource_id": "SMIC_TOKYO_GENOME_2025",
                "target_seqid": "JBLTJD010000033.1",
                "target_start": "100",
                "target_end": "300",
                "strand": "+",
                "feature_id": "MP1",
                "feature_type": "mRNA",
                "query_start": "1",
                "query_end": "100",
                "query_length_aa": "100",
                "query_coverage": "1",
                "cds_feature_count": "4",
                "exon_feature_count": "0",
                "frameshift_or_stop_flag": "False",
                "identity": "0.99",
                "positive": "1",
                "score": "100",
                "rank": "1",
                "alignment_status": "HIGH_COVERAGE_NO_DISRUPTION",
                "raw_attributes": "ID=MP1",
                "notes": "test",
            }
        ],
        PHASE4B_ALIGNMENT_HIT_COLUMNS,
    )
    raw_gff.write_text(
        "##gff-version 3\n"
        "##PAF\tFTH1B__q1\t100\t0\t100\t+\tJBLTJD010000033.1\t1000\t100\t300\t300\t300\t0\tfs:i:0\tst:i:0\n"
        "JBLTJD010000033.1\tminiprot\tmRNA\t100\t300\t99\t+\t.\tID=MP1;Target=FTH1B__q1 1 100\n",
        encoding="utf-8",
    )
    write_tsv(
        domains,
        [
            {
                "mechanism": "Ferroptosis, iron handling, and oxidative stress",
                "gene_symbol": "FTH1B",
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "gene_id": "prot1",
                "orthogroup_id": "NOT_ASSESSED",
                "representative_protein_id": "prot1",
                "representative_fasta_id": "FTH1B__q1",
                "representative_length_aa": "100",
                "domain_validation_status": "DOMAIN_SUPPORTED",
                "best_pfam_accession": "PF00210",
                "best_pfam_name": "Ferritin",
            }
        ],
        DOMAIN_INTEGRITY_COLUMNS,
    )
    write_tsv(
        target_regions,
        [
            {
                "gene_symbol": "FTH1B",
                "source_coordinate": "scaffold_33",
                "requested_seqid": "scaffold_33",
                "target_fasta_header": "JBLTJD010000033.1 test scaffold_33",
                "target_fasta_path": "target.fna",
                "extraction_status": "EXTRACTED",
                "notes": "test",
            }
        ],
        PHASE4B_TARGET_REGION_COLUMNS,
    )

    locus_rows, gene_rows = review_loci(hits, raw_gff, domains, target_regions, config)
    assert locus_rows[0]["domain_validation_status"] == "DOMAIN_SUPPORTED"
    assert locus_rows[0]["miniprot_disruption_status"] == "NO_PARSED_FS_OR_ST_TAG"
    assert gene_rows[0]["reviewed_locus_count"] == "1"


def test_tp53_no_alignment_summary_is_annotation_uncertainty():
    row = tp53_summary_row(None, [], [], ready=True, full_coverage=0.70, partial_coverage=0.20)
    assert row["genome_search_status"] == "NO_TP53_ALIGNMENT_UNDER_CURRENT_QUERY"
    assert row["candidate_classification"] == "annotation_uncertainty"
    assert "absence" in row["claim_language_guardrail"]


def test_tp53_disrupted_high_coverage_summary_is_not_generic_partial():
    row = tp53_summary_row(
        None,
        [
            {
                "target_seqid": "scaffold_247",
                "target_start": "10",
                "target_end": "100",
                "query_coverage": "0.96",
                "frameshift_or_stop_flag": "True",
            }
        ],
        [],
        ready=True,
        full_coverage=0.70,
        partial_coverage=0.20,
    )
    assert row["genome_search_status"] == "POSSIBLE_DISRUPTED_TP53_ALIGNMENT"
    assert row["candidate_classification"] == "possible_disrupted_or_divergent_candidate_locus"
    assert row["candidate_loci"] == "scaffold_247:10-100"


def test_tp53_forward_target_scaffolds_uses_coordinate_seqids(tmp_path: Path):
    reciprocal = tmp_path / "reciprocal.tsv"
    write_tsv(
        reciprocal,
        [
            {
                "gene_symbol": "TP53",
                "target_protein_id": "p53like",
                "target_description": "cellular tumor antigen p53-like",
                "target_gene_id": "gene-gs_1",
                "target_coordinate_summary": "gene-gs_1:scaffold_247:10-100:-",
                "source_species_id": "rtyp",
                "source_resource_id": "RTYP_REFSEQ_2022",
                "query_protein_id": "XP_1",
                "reciprocal_top_protein_id": "other",
                "reciprocal_top_description": "other",
                "reciprocal_top_evalue": "1e-5",
                "reciprocal_top_bitscore": "50",
                "reciprocal_gene_pattern_match": "false",
                "reciprocal_selected_query_match": "false",
                "reciprocal_support_status": "NOT_RECIPROCAL",
                "notes": "test",
            }
        ],
        PHASE4_RECIPROCAL_HIT_COLUMNS,
    )
    scaffolds, source_rows = tp53_forward_target_scaffolds(reciprocal)
    assert scaffolds == {"TP53": {"scaffold_247"}}
    assert source_rows[0]["requested_seqid"] == "scaffold_247"
