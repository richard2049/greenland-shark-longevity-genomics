from pathlib import Path

import yaml

from greenland_shark_longevity.phase7c_expression_reference import (
    parse_gff_models,
    run_phase7c,
    scaffold_aliases,
)
from greenland_shark_longevity.schemas import (
    PHASE7B_CANDIDATE_QUANT_MAP_COLUMNS,
    PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS,
    PHASE7C_REFERENCE_QC_COLUMNS,
    PHASE7C_TX2GENE_COLUMNS,
)
from greenland_shark_longevity.utils import read_tsv, write_tsv


def write_phase7c_fixture(tmp_path: Path) -> Path:
    gff = tmp_path / "annotation.gff"
    genome = tmp_path / "genome.fna"
    candidate_map = tmp_path / "phase7b_map.tsv"
    config = tmp_path / "config.yaml"

    gff.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "scaffold_1\tGnomon\tgene\t1\t12\t.\t+\t.\tID=gene-gs_000001;gene=ERCC1;locus_tag=gs_000001",
                "scaffold_1\tGnomon\tmRNA\t1\t12\t.\t+\t.\tID=rna-gs_000001-R1;Parent=gene-gs_000001;gene=ERCC1;locus_tag=gs_000001;product=ERCC1 transcript",
                "scaffold_1\tGnomon\texon\t1\t3\t.\t+\t.\tID=exon1;Parent=rna-gs_000001-R1",
                "scaffold_1\tGnomon\texon\t10\t12\t.\t+\t.\tID=exon2;Parent=rna-gs_000001-R1",
                "scaffold_2\tGnomon\tgene\t2\t4\t.\t-\t.\tID=gene-gs_000002;gene=FTH1B;locus_tag=gs_000002",
                "scaffold_2\tGnomon\tmRNA\t2\t4\t.\t-\t.\tID=rna-gs_000002-R1;Parent=gene-gs_000002;gene=FTH1B;locus_tag=gs_000002;product=ferritin transcript",
                "scaffold_2\tGnomon\texon\t2\t4\t.\t-\t.\tID=exon3;Parent=rna-gs_000002-R1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome.write_text(
        ">JBLTJD000001.1 test scaffold_1 whole genome shotgun sequence\nAAACCCGGGTTT\n"
        ">JBLTJD000002.1 test scaffold_2 whole genome shotgun sequence\nAACCGG\n",
        encoding="utf-8",
    )
    write_tsv(
        candidate_map,
        [
            {
                "mechanism": "DNA repair/genome stability",
                "gene_symbol": "ERCC1",
                "quantification_gene_ids": "gene-gs_000001",
                "quantification_strategy_status": "READY_FOR_FIRST_PASS_CANDIDATE_QUANTIFICATION",
            },
            {
                "mechanism": "Ferroptosis, iron handling, and oxidative stress",
                "gene_symbol": "FTH1B",
                "quantification_gene_ids": "gene-gs_000002",
                "quantification_strategy_status": "TARGETED_LOCUS_VALIDATION_REQUIRED_AFTER_FIRST_PASS_QUANTIFICATION",
            },
            {
                "mechanism": "p53 pathway",
                "gene_symbol": "TP53",
                "quantification_gene_ids": "NOT_ASSESSED",
                "quantification_strategy_status": "NOT_READY_NO_CANDIDATE_GENE_MODEL_MAP",
            },
        ],
        PHASE7B_CANDIDATE_QUANT_MAP_COLUMNS,
    )
    config.write_text(
        yaml.safe_dump(
            {
                "phase7c_expression_reference_validation": {
                    "annotation_gff": str(gff),
                    "genome_fasta": str(genome),
                    "phase7b_candidate_quant_map": str(candidate_map),
                    "transcript_fasta_output": str(tmp_path / "transcripts.fna"),
                    "tx2gene_output": str(tmp_path / "tx2gene.tsv"),
                    "reference_qc_output": str(tmp_path / "reference_qc.tsv"),
                    "candidate_reference_validation_output": str(tmp_path / "candidate_validation.tsv"),
                    "report_output": str(tmp_path / "report.md"),
                    "transcript_feature_types": ["mRNA"],
                    "interval_feature_preference": "exon",
                    "fallback_interval_feature": "CDS",
                }
            }
        ),
        encoding="utf-8",
    )
    return config


def test_scaffold_aliases_detects_accession_and_scaffold_name():
    aliases = scaffold_aliases(">JBLTJD000001.1 Somniosus scaffold_1, whole genome shotgun sequence")

    assert "JBLTJD000001.1" in aliases
    assert "scaffold_1" in aliases


def test_parse_gff_models_uses_exon_intervals(tmp_path: Path):
    config_path = write_phase7c_fixture(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    phase7c = config["phase7c_expression_reference_validation"]

    models = parse_gff_models(Path(phase7c["annotation_gff"]), {"mRNA"}, "exon", "CDS")

    assert models["rna-gs_000001-R1"].gene_id == "gene-gs_000001"
    assert models["rna-gs_000001-R1"].intervals == [(1, 3), (10, 12)]
    assert models["rna-gs_000002-R1"].strand == "-"


def test_phase7c_outputs_reference_and_candidate_validation(tmp_path: Path):
    config_path = write_phase7c_fixture(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    phase7c = config["phase7c_expression_reference_validation"]

    run_phase7c(
        config_path,
        Path(phase7c["transcript_fasta_output"]),
        Path(phase7c["tx2gene_output"]),
        Path(phase7c["reference_qc_output"]),
        Path(phase7c["candidate_reference_validation_output"]),
        Path(phase7c["report_output"]),
    )

    fasta_text = Path(phase7c["transcript_fasta_output"]).read_text(encoding="utf-8")
    assert "AAATTT" in fasta_text
    assert "GGT" in fasta_text

    tx2gene = read_tsv(Path(phase7c["tx2gene_output"]), PHASE7C_TX2GENE_COLUMNS)
    assert len(tx2gene) == 2
    assert tx2gene[0]["reference_status"] == "TRANSCRIPT_REFERENCE_CONSTRUCTED"

    qc = read_tsv(Path(phase7c["reference_qc_output"]), PHASE7C_REFERENCE_QC_COLUMNS)[0]
    assert qc["transcript_count"] == "2"
    assert qc["validation_status"] == "REFERENCE_CONSTRUCTED_AND_INPUTS_VALIDATED"

    validation = read_tsv(
        Path(phase7c["candidate_reference_validation_output"]),
        PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS,
    )
    by_gene = {row["gene_symbol"]: row for row in validation}
    assert by_gene["ERCC1"]["phase7c_quantification_readiness"] == "REFERENCE_READY_FOR_FIRST_PASS_QUANTIFICATION"
    assert by_gene["FTH1B"]["phase7c_quantification_readiness"] == "REFERENCE_PRESENT_TARGETED_LOCUS_VALIDATION_REQUIRED"
    assert by_gene["TP53"]["phase7c_quantification_readiness"] == "NOT_QUANTIFIABLE_NO_GENE_MODEL_MAP"
    assert "does not download reads" in Path(phase7c["report_output"]).read_text(encoding="utf-8")
