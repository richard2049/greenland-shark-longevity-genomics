from pathlib import Path

from greenland_shark_longevity.genome_validation import (
    QueryRecord,
    alignment_rows,
    build_query_records,
    extract_target_scaffolds,
    parse_miniprot_gff,
    summary_rows,
)
from greenland_shark_longevity.schemas import PHASE4_RESCUE_QUERY_COLUMNS, PHASE4_RESCUE_SUMMARY_COLUMNS
from greenland_shark_longevity.utils import write_tsv


def test_phase4b_queries_include_focal_rescue_and_unresolved_comparator(tmp_path: Path):
    summary = tmp_path / "summary.tsv"
    query_inventory = tmp_path / "queries.tsv"
    rescue_fasta = tmp_path / "rescue.faa"
    comparator_fasta = tmp_path / "comparator.faa"
    write_tsv(
        summary,
        [
            {
                "gene_symbol": "H1F0",
                "mechanism": "Chromatin and epigenome-related regulation",
                "protein_rescue_status": "PROTEIN_LEVEL_RESCUE_CANDIDATE",
                "candidate_protein_ids": "protA",
            },
            {
                "gene_symbol": "TP53",
                "mechanism": "p53 pathway",
                "protein_rescue_status": "NO_RECIPROCAL_PROTEIN_RESCUE",
                "candidate_protein_ids": "NOT_ASSESSED",
            },
        ],
        PHASE4_RESCUE_SUMMARY_COLUMNS,
    )
    write_tsv(
        query_inventory,
        [
            {
                "gene_symbol": "TP53",
                "mechanism": "p53 pathway",
                "query_species_id": "scan",
                "query_resource_id": "SCAN_REFSEQ_2026",
                "query_protein_id": "qTP53",
                "query_description": "cellular tumor antigen p53 [test]",
                "query_length_aa": "4",
                "source_fasta": str(comparator_fasta),
                "selection_status": "SELECTED_FOR_RESCUE_QUERY",
            }
        ],
        PHASE4_RESCUE_QUERY_COLUMNS,
    )
    rescue_fasta.write_text(">protA histone H1.0-B-like\nMAAA\n", encoding="utf-8")
    comparator_fasta.write_text(">qTP53 cellular tumor antigen p53 [test]\nMBBB\n", encoding="utf-8")

    queries = build_query_records(summary, query_inventory, rescue_fasta)
    assert {query.gene_symbol for query in queries} == {"H1F0", "TP53"}
    assert {query.query_source_type for query in queries} == {
        "focal_rescue_candidate",
        "comparator_query_for_unresolved_gene",
    }


def test_parse_miniprot_gff_and_classify_high_coverage(tmp_path: Path):
    gff = tmp_path / "aln.gff"
    query = QueryRecord(
        gene_symbol="RAD51",
        mechanism="DNA repair/genome stability",
        query_id="RAD51__focal_rescue_candidate__prot1",
        query_source_type="focal_rescue_candidate",
        original_protein_id="prot1",
        original_description="DNA repair protein RAD51 homolog 1",
        source_species_id="smic",
        source_resource_id="SMIC_TOKYO_GENOME_2025",
        source_fasta="rescue.faa",
        sequence="M" * 100,
        selection_reason="test",
    )
    gff.write_text(
        "##gff-version 3\n"
        "scaffold_1\tminiprot\tmRNA\t10\t300\t99\t+\t.\tID=aln1;Target=RAD51__focal_rescue_candidate__prot1 1 90;Identity=0.95;Rank=1\n"
        "scaffold_1\tminiprot\tCDS\t10\t100\t.\t+\t0\tParent=aln1\n",
        encoding="utf-8",
    )
    hits = parse_miniprot_gff(gff, {query.query_id: query})
    rows = alignment_rows(hits, full_coverage=0.70, partial_coverage=0.20)
    assert len(hits) == 1
    assert hits[0].query_coverage == 0.9
    assert rows[0]["alignment_status"] == "HIGH_COVERAGE_NO_DISRUPTION"


def test_parse_miniprot_terminal_stop_codon_is_not_disruption(tmp_path: Path):
    gff = tmp_path / "aln.gff"
    query = QueryRecord(
        gene_symbol="FTH1B",
        mechanism="Ferroptosis, iron handling, and oxidative stress",
        query_id="FTH1B__focal_rescue_candidate__prot1",
        query_source_type="focal_rescue_candidate",
        original_protein_id="prot1",
        original_description="ferritin heavy chain",
        source_species_id="smic",
        source_resource_id="SMIC_TOKYO_GENOME_2025",
        source_fasta="rescue.faa",
        sequence="M" * 100,
        selection_reason="test",
    )
    gff.write_text(
        "##gff-version 3\n"
        "##PAF\tFTH1B__focal_rescue_candidate__prot1\t100\t0\t100\t+\tscaffold_1\t1000\t10\t310\t300\t300\t0\tfs:i:0\tst:i:0\n"
        "scaffold_1\tminiprot\tmRNA\t10\t310\t99\t+\t.\tID=aln1;Target=FTH1B__focal_rescue_candidate__prot1 1 100;Identity=0.95;Rank=1\n"
        "scaffold_1\tminiprot\tCDS\t10\t307\t.\t+\t0\tParent=aln1\n"
        "scaffold_1\tminiprot\tstop_codon\t308\t310\t.\t+\t0\tParent=aln1\n",
        encoding="utf-8",
    )
    hits = parse_miniprot_gff(gff, {query.query_id: query})
    rows = alignment_rows(hits, full_coverage=0.70, partial_coverage=0.20)
    assert len(hits) == 1
    assert hits[0].frameshift_or_stop_flag is False
    assert rows[0]["alignment_status"] == "HIGH_COVERAGE_NO_DISRUPTION"


def test_parse_miniprot_paf_frameshift_or_stop_tag_is_disruption(tmp_path: Path):
    gff = tmp_path / "aln.gff"
    query = QueryRecord(
        gene_symbol="FTH1B",
        mechanism="Ferroptosis, iron handling, and oxidative stress",
        query_id="FTH1B__focal_rescue_candidate__prot1",
        query_source_type="focal_rescue_candidate",
        original_protein_id="prot1",
        original_description="ferritin heavy chain",
        source_species_id="smic",
        source_resource_id="SMIC_TOKYO_GENOME_2025",
        source_fasta="rescue.faa",
        sequence="M" * 100,
        selection_reason="test",
    )
    gff.write_text(
        "##gff-version 3\n"
        "##PAF\tFTH1B__focal_rescue_candidate__prot1\t100\t0\t100\t+\tscaffold_1\t1000\t10\t310\t300\t300\t0\tfs:i:1\tst:i:0\n"
        "scaffold_1\tminiprot\tmRNA\t10\t310\t99\t+\t.\tID=aln1;Target=FTH1B__focal_rescue_candidate__prot1 1 100;Identity=0.95;Rank=1\n",
        encoding="utf-8",
    )
    hits = parse_miniprot_gff(gff, {query.query_id: query})
    rows = alignment_rows(hits, full_coverage=0.70, partial_coverage=0.20)
    assert len(hits) == 1
    assert hits[0].frameshift_or_stop_flag is True
    assert rows[0]["alignment_status"] == "POSSIBLE_DISRUPTION"


def test_phase4b_summary_missing_inputs_is_not_absence(tmp_path: Path):
    summary = tmp_path / "summary.tsv"
    write_tsv(
        summary,
        [
            {
                "gene_symbol": "TP53",
                "mechanism": "p53 pathway",
                "protein_rescue_status": "NO_RECIPROCAL_PROTEIN_RESCUE",
            }
        ],
        PHASE4_RESCUE_SUMMARY_COLUMNS,
    )
    rows = summary_rows(summary, queries=[], hits=[], ready=False, full_coverage=0.70, partial_coverage=0.20)
    assert rows[0]["genome_validation_status"] == "NOT_RUN_REQUIRES_INPUTS"
    assert rows[0]["candidate_classification"] == "not_assessable"
    assert "absent" in rows[0]["claim_language_guardrail"]


def test_extract_target_scaffolds_matches_exact_scaffold_alias(tmp_path: Path):
    genome = tmp_path / "genome.fa"
    target = tmp_path / "target.fa"
    genome.write_text(
        ">ACC12 species scaffold_12, whole genome shotgun sequence\nAAAA\n"
        ">ACC120 species scaffold_120, whole genome shotgun sequence\nCCCC\n"
        ">ACC33 species scaffold_33, whole genome shotgun sequence\nGGGG\n",
        encoding="utf-8",
    )
    found, rows = extract_target_scaffolds(
        genome,
        target,
        {"RAD51": {"scaffold_12"}, "H1F0": {"scaffold_33"}},
    )
    assert found == {"scaffold_12", "scaffold_33"}
    text = target.read_text(encoding="utf-8")
    assert "scaffold_12," in text
    assert "scaffold_33," in text
    assert "scaffold_120" not in text
    assert {row["extraction_status"] for row in rows} == {"EXTRACTED"}
