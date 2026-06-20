from pathlib import Path

import yaml

from greenland_shark_longevity.phase7b_quantification_strategy import (
    build_candidate_quant_map,
    build_reference_strategy_rows,
    parse_cds_headers,
    write_phase7b_outputs,
)
from greenland_shark_longevity.schemas import (
    CANDIDATE_ISOFORM_AUDIT_COLUMNS,
    COPY_NUMBER_COLUMNS,
    PHASE4E_GENE_HARDENING_COLUMNS,
    PHASE4E_LOCUS_HARDENING_COLUMNS,
    PHASE7_CANDIDATE_EXPRESSION_PLAN_COLUMNS,
    PHASE7_RNASEQ_READINESS_COLUMNS,
    PHASE7B_CANDIDATE_QUANT_MAP_COLUMNS,
    PHASE7B_REFERENCE_STRATEGY_COLUMNS,
    PHASE7B_RUN_PLAN_COLUMNS,
)
from greenland_shark_longevity.utils import read_tsv, write_tsv


def write_phase7b_fixture(tmp_path: Path) -> Path:
    plan = tmp_path / "phase7a_plan.tsv"
    readiness = tmp_path / "readiness.tsv"
    copy = tmp_path / "copy.tsv"
    isoforms = tmp_path / "isoforms.tsv"
    phase4e_gene = tmp_path / "phase4e_gene.tsv"
    phase4e_locus = tmp_path / "phase4e_locus.tsv"
    cds = tmp_path / "cds.fna"
    gff = tmp_path / "annotation.gff"
    genome = tmp_path / "genome.fna"
    config = tmp_path / "config.yaml"

    write_tsv(
        plan,
        [
            {
                "mechanism": "DNA repair/genome stability",
                "gene_symbol": "ERCC1",
                "usable_rnaseq_run_count": "3",
                "usable_tissues": "retina",
                "current_expression_status": "NOT_QUANTIFIED_PHASE7A_METADATA_ONLY",
            },
            {
                "mechanism": "Ferroptosis, iron handling, and oxidative stress",
                "gene_symbol": "FTH1B",
                "usable_rnaseq_run_count": "3",
                "usable_tissues": "retina",
                "current_expression_status": "NOT_QUANTIFIED_PHASE7A_METADATA_ONLY",
            },
            {
                "mechanism": "p53 pathway",
                "gene_symbol": "TP53",
                "usable_rnaseq_run_count": "3",
                "usable_tissues": "retina",
                "current_expression_status": "NOT_QUANTIFIED_PHASE7A_METADATA_ONLY",
            },
        ],
        PHASE7_CANDIDATE_EXPRESSION_PLAN_COLUMNS,
    )
    write_tsv(
        readiness,
        [
            {
                "resource_id": "SMIC_RETINA",
                "run": "SRR1",
                "included_for_candidate_expression_plan": "true",
                "expression_readiness_status": "RNA_SEQ_METADATA_READY_FOR_CANDIDATE_EXPRESSION_AUDIT",
            }
        ],
        PHASE7_RNASEQ_READINESS_COLUMNS,
    )
    write_tsv(
        copy,
        [
            {
                "mechanism": "DNA repair/genome stability",
                "gene_symbol": "ERCC1",
                "gene_ids": "gene-gs_000001",
                "copy_count": "1",
                "mapping_status": "ANNOTATION_SYMBOL_MATCH_ORTHOGROUP_MAPPED",
            },
            {
                "mechanism": "Ferroptosis, iron handling, and oxidative stress",
                "gene_symbol": "FTH1B",
                "gene_ids": "NOT_ASSESSED",
                "copy_count": "0",
                "mapping_status": "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH",
            },
            {
                "mechanism": "p53 pathway",
                "gene_symbol": "TP53",
                "gene_ids": "NOT_ASSESSED",
                "copy_count": "0",
                "mapping_status": "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH",
            },
        ],
        COPY_NUMBER_COLUMNS,
    )
    write_tsv(
        isoforms,
        [
            {
                "mechanism": "DNA repair/genome stability",
                "gene_symbol": "ERCC1",
                "gene_id": "gene-gs_000001",
                "isoform_count": "1",
                "representative_protein_id": "gnl|WGS_ZZZZ|gs_000001-P1",
            }
        ],
        CANDIDATE_ISOFORM_AUDIT_COLUMNS,
    )
    write_tsv(
        phase4e_gene,
        [
            {
                "gene_symbol": "FTH1B",
                "mechanism": "Ferroptosis, iron handling, and oxidative stress",
                "phase4e_hardened_status": "HARDENED_CANDIDATE_CLUSTER_REQUIRES_FERRITIN_FAMILY_RESOLUTION",
                "evidence_tier_recommendation": "Tier 2",
                "artifact_risk": "high",
            },
            {
                "gene_symbol": "TP53",
                "mechanism": "p53 pathway",
                "phase4e_hardened_status": "UNRESOLVED_P53_FAMILY_CANDIDATE_REQUIRES_ADDITIONAL_PHASE4_VALIDATION",
                "evidence_tier_recommendation": "Artifact/uncertain",
                "artifact_risk": "high",
            },
        ],
        PHASE4E_GENE_HARDENING_COLUMNS,
    )
    write_tsv(
        phase4e_locus,
        [
            {
                "gene_symbol": "FTH1B",
                "mechanism": "Ferroptosis, iron handling, and oxidative stress",
                "locus_cluster_id": "FTH1B_L001",
                "overlapping_gene_ids": "gene-gs_000002",
            },
            {
                "gene_symbol": "TP53",
                "mechanism": "p53 pathway",
                "locus_cluster_id": "TP53_L001",
                "overlapping_gene_ids": "gene-gs_000003",
            },
        ],
        PHASE4E_LOCUS_HARDENING_COLUMNS,
    )
    cds.write_text(
        ">lcl|scaffold_1_cds_gs_000001-P1_1 [gene=ERCC1] [locus_tag=gs_000001] [protein_id=gs_000001-P1]\nATGC\n"
        ">lcl|scaffold_2_cds_gs_000002-P1_1 [gene=FTH1B] [locus_tag=gs_000002] [protein_id=gs_000002-P1]\nATGC\n"
        ">lcl|scaffold_3_cds_gs_000003-P1_1 [gene=TP53] [locus_tag=gs_000003] [protein_id=gs_000003-P1]\nATGC\n",
        encoding="utf-8",
    )
    gff.write_text("##gff-version 3\n", encoding="utf-8")
    genome.write_text(">scaffold_1\nATGC\n", encoding="utf-8")
    config.write_text(
        yaml.safe_dump(
            {
                "phase7_rnaseq": {
                    "registered_runs": [
                        {"run": "SRR1", "library_strategy": "RNA-Seq"},
                        {"run": "SRR2", "library_strategy": "WGS"},
                    ]
                },
                "phase7b_expression_quantification_strategy": {
                    "phase7a_candidate_expression_plan": str(plan),
                    "phase7a_readiness": str(readiness),
                    "candidate_copy_number": str(copy),
                    "candidate_isoform_audit": str(isoforms),
                    "phase4e_gene_summary": str(phase4e_gene),
                    "phase4e_locus_review": str(phase4e_locus),
                    "annotation_gff": str(gff),
                    "genome_fasta": str(genome),
                    "cds_fasta": str(cds),
                    "transcriptome_fasta_planned": str(tmp_path / "transcripts.fna"),
                    "tx2gene_planned": str(tmp_path / "tx2gene.tsv"),
                    "salmon_index_planned": str(tmp_path / "salmon_index"),
                    "candidate_quantification_matrix_planned": str(tmp_path / "matrix.tsv"),
                    "high_ambiguity_genes": ["FTH1B", "TP53"],
                },
            }
        ),
        encoding="utf-8",
    )
    return config


def test_parse_cds_headers_extracts_locus_and_protein(tmp_path: Path):
    cds = tmp_path / "cds.fna"
    cds.write_text(
        ">seq1 [gene=ERCC1] [locus_tag=gs_000001] [protein_id=gs_000001-P1]\nATGC\n",
        encoding="utf-8",
    )

    rows = parse_cds_headers(cds)

    assert rows[0]["gene_symbol"] == "ERCC1"
    assert rows[0]["locus_tag"] == "gs_000001"
    assert rows[0]["protein_id"] == "gs_000001-P1"


def test_candidate_quantification_map_classifies_ready_and_high_risk(tmp_path: Path):
    config_path = write_phase7b_fixture(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    rows = build_candidate_quant_map(config)
    by_gene = {row["gene_symbol"]: row for row in rows}

    assert by_gene["ERCC1"]["quantification_strategy_status"] == "READY_FOR_FIRST_PASS_CANDIDATE_QUANTIFICATION"
    assert by_gene["ERCC1"]["cds_reference_sequence_ids"] == "lcl|scaffold_1_cds_gs_000001-P1_1"
    assert by_gene["FTH1B"]["quantification_strategy_status"] == "TARGETED_LOCUS_VALIDATION_REQUIRED_AFTER_FIRST_PASS_QUANTIFICATION"
    assert by_gene["FTH1B"]["phase4e_candidate_gene_ids"] == "gene-gs_000002"
    assert by_gene["TP53"]["quantification_strategy_status"] == "NOT_READY_TARGETED_LOCUS_VALIDATION_REQUIRED_BEFORE_EXPRESSION_INTERPRETATION"


def test_reference_strategy_and_cli_outputs(tmp_path: Path):
    config_path = write_phase7b_fixture(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    ref_out = tmp_path / "reference.tsv"
    map_out = tmp_path / "map.tsv"
    plan_out = tmp_path / "run_plan.tsv"
    report = tmp_path / "report.md"

    reference_rows = build_reference_strategy_rows(config)
    write_phase7b_outputs(config_path, ref_out, map_out, plan_out, report)

    assert reference_rows[0]["strategy_id"] == "PHASE7B_PRIMARY_SALMON_SELECTIVE_ALIGNMENT"
    assert len(read_tsv(ref_out, PHASE7B_REFERENCE_STRATEGY_COLUMNS)) == 3
    assert len(read_tsv(map_out, PHASE7B_CANDIDATE_QUANT_MAP_COLUMNS)) == 3
    assert read_tsv(plan_out, PHASE7B_RUN_PLAN_COLUMNS)
    assert "does not download reads" in report.read_text(encoding="utf-8")
