from pathlib import Path

import yaml

from greenland_shark_longevity.phase7d_rnaseq_quantification import (
    build_candidate_expression_matrix,
    build_fastq_qc_rows,
    build_raw_read_intake_rows,
    run_salmon_quantification,
    run_phase7d,
)
from greenland_shark_longevity.schemas import (
    PHASE7_RNASEQ_MANIFEST_COLUMNS,
    PHASE7_RNASEQ_READINESS_COLUMNS,
    PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS,
    PHASE7C_TX2GENE_COLUMNS,
    PHASE7D_CANDIDATE_EXPRESSION_MATRIX_COLUMNS,
    PHASE7D_FASTQ_QC_COLUMNS,
    PHASE7D_RAW_READ_INTAKE_COLUMNS,
    PHASE7D_SALMON_PREFLIGHT_COLUMNS,
    PHASE7D_SALMON_QUANT_SUMMARY_COLUMNS,
)
from greenland_shark_longevity.utils import read_tsv, write_tsv


def write_phase7d_fixture(tmp_path: Path, with_fastq: bool = False) -> Path:
    raw_dir = tmp_path / "raw_reads"
    raw_dir.mkdir()
    transcript_fasta = tmp_path / "transcripts.fna"
    tx2gene = tmp_path / "tx2gene.tsv"
    manifest = tmp_path / "rnaseq_manifest.tsv"
    readiness = tmp_path / "readiness.tsv"
    candidate_validation = tmp_path / "candidate_reference_validation.tsv"
    config = tmp_path / "config.yaml"

    transcript_fasta.write_text(">tx1 gene=RAD51\nATGC\n>tx2 gene=RAD51\nATGC\n>tx3 gene=TP53\nATGC\n", encoding="utf-8")
    write_tsv(
        tx2gene,
        [
            {
                "transcript_id": "tx1",
                "gene_id": "gene-rad51",
                "gene_symbol": "RAD51",
                "locus_tag": "rad51",
                "scaffold": "scaffold_1",
                "strand": "+",
                "feature_type": "exon",
                "exon_count": "1",
                "transcript_length_bp": "4",
                "reference_status": "TRANSCRIPT_REFERENCE_CONSTRUCTED",
            },
            {
                "transcript_id": "tx2",
                "gene_id": "gene-rad51",
                "gene_symbol": "RAD51",
                "locus_tag": "rad51",
                "scaffold": "scaffold_1",
                "strand": "+",
                "feature_type": "exon",
                "exon_count": "1",
                "transcript_length_bp": "4",
                "reference_status": "TRANSCRIPT_REFERENCE_CONSTRUCTED",
            },
            {
                "transcript_id": "tx3",
                "gene_id": "gene-tp53",
                "gene_symbol": "TP53",
                "locus_tag": "tp53",
                "scaffold": "scaffold_2",
                "strand": "-",
                "feature_type": "exon",
                "exon_count": "1",
                "transcript_length_bp": "4",
                "reference_status": "TRANSCRIPT_REFERENCE_CONSTRUCTED",
            },
        ],
        PHASE7C_TX2GENE_COLUMNS,
    )
    write_tsv(
        manifest,
        [
            {
                "resource_id": "SMIC_RETINA_PRJNA1246101_2026",
                "bioproject": "PRJNA1246101",
                "sra_study": "SRP576249",
                "experiment": "SRX1",
                "run": "SRR1",
                "biosample": "SAMN1",
                "sample": "SRS1",
                "sample_name": "Greenland_shark_retina_1",
                "scientific_name": "Somniosus microcephalus",
                "taxon_id": "191813",
                "sex": "female",
                "tissue": "retina",
                "library_strategy": "RNA-Seq",
                "library_source": "TRANSCRIPTOMIC",
                "library_selection": "PolyA",
                "library_layout": "PAIRED",
                "platform": "ILLUMINA",
                "model": "NovaSeq",
                "spots": "10",
                "bases": "2000",
                "avg_length": "100",
                "size_mb": "1",
                "release_date": "2026-01-01",
                "load_date": "2026-01-02",
                "download_path": "NOT_ASSESSED",
                "source_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1246101",
                "retrieval_date": "2026-06-12",
                "local_path": str(raw_dir / "SRR1.sra"),
                "download_status": "RAW_READS_NOT_DOWNLOADED_PHASE7A_METADATA_ONLY",
                "metadata_status": "RUNINFO_METADATA_REGISTERED",
                "usage_notes": "fixture",
            },
            {
                "resource_id": "SMIC_RETINA_PRJNA1246101_2026",
                "bioproject": "PRJNA1246101",
                "sra_study": "SRP576249",
                "experiment": "SRX2",
                "run": "SRR2",
                "biosample": "SAMN2",
                "sample": "SRS2",
                "sample_name": "Greenland_shark_wgs",
                "scientific_name": "Somniosus microcephalus",
                "taxon_id": "191813",
                "sex": "female",
                "tissue": "retina",
                "library_strategy": "WGS",
                "library_source": "GENOMIC",
                "library_selection": "other",
                "library_layout": "PAIRED",
                "platform": "ILLUMINA",
                "model": "NovaSeq",
                "spots": "10",
                "bases": "2000",
                "avg_length": "100",
                "size_mb": "1",
                "release_date": "2026-01-01",
                "load_date": "2026-01-02",
                "download_path": "NOT_ASSESSED",
                "source_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1246101",
                "retrieval_date": "2026-06-12",
                "local_path": str(raw_dir / "SRR2.sra"),
                "download_status": "RAW_READS_NOT_DOWNLOADED_PHASE7A_METADATA_ONLY",
                "metadata_status": "RUNINFO_METADATA_REGISTERED",
                "usage_notes": "fixture",
            },
        ],
        PHASE7_RNASEQ_MANIFEST_COLUMNS,
    )
    write_tsv(
        readiness,
        [
            {
                "resource_id": "SMIC_RETINA_PRJNA1246101_2026",
                "run": "SRR1",
                "experiment": "SRX1",
                "biosample": "SAMN1",
                "sample_name": "Greenland_shark_retina_1",
                "tissue": "retina",
                "library_strategy": "RNA-Seq",
                "library_source": "TRANSCRIPTOMIC",
                "library_selection": "PolyA",
                "library_layout": "PAIRED",
                "included_for_candidate_expression_plan": "true",
            },
            {
                "resource_id": "SMIC_RETINA_PRJNA1246101_2026",
                "run": "SRR2",
                "experiment": "SRX2",
                "biosample": "SAMN2",
                "sample_name": "Greenland_shark_wgs",
                "tissue": "retina",
                "library_strategy": "WGS",
                "library_source": "GENOMIC",
                "library_selection": "other",
                "library_layout": "PAIRED",
                "included_for_candidate_expression_plan": "false",
            },
        ],
        PHASE7_RNASEQ_READINESS_COLUMNS,
    )
    write_tsv(
        candidate_validation,
        [
            {
                "mechanism": "DNA repair/genome stability",
                "gene_symbol": "RAD51",
                "matched_transcript_ids": "tx1;tx2",
                "matched_transcript_count": "2",
                "phase7c_quantification_readiness": "REFERENCE_READY_FOR_FIRST_PASS_QUANTIFICATION",
            },
            {
                "mechanism": "p53 pathway",
                "gene_symbol": "TP53",
                "matched_transcript_ids": "tx3",
                "matched_transcript_count": "1",
                "phase7c_quantification_readiness": "REFERENCE_PRESENT_TARGETED_LOCUS_VALIDATION_REQUIRED",
            },
        ],
        PHASE7C_CANDIDATE_REFERENCE_VALIDATION_COLUMNS,
    )
    if with_fastq:
        (raw_dir / "SRR1_1.fastq").write_text("@r1\nACGT\n+\nIIII\n@r2\nGGNN\n+\nHHHH\n", encoding="utf-8")
        (raw_dir / "SRR1_2.fastq").write_text("@r1\nTTAA\n+\nIIII\n@r2\nCCCC\n+\nHHHH\n", encoding="utf-8")
    config.write_text(
        yaml.safe_dump(
            {
                "phase7d_rnaseq_quantification": {
                    "phase7a_manifest": str(manifest),
                    "phase7a_readiness": str(readiness),
                    "phase7c_tx2gene": str(tx2gene),
                    "phase7c_candidate_reference_validation": str(candidate_validation),
                    "transcript_fasta": str(transcript_fasta),
                    "raw_read_local_dir": str(raw_dir),
                    "salmon_executable": "salmon_missing_for_test",
                    "salmon_index_dir": str(tmp_path / "salmon_index"),
                    "salmon_quant_dir": str(tmp_path / "salmon_quant"),
                    "logs_dir": str(tmp_path / "logs"),
                    "expected_fastq_patterns": ["{run}_1.fastq", "{run}_2.fastq", "{run}_R1.fastq.gz", "{run}_R2.fastq.gz"],
                    "fastq_qc_max_records": 10,
                    "build_salmon_index_if_missing": False,
                    "run_salmon_if_ready": True,
                    "reuse_existing_salmon_quant_if_present": True,
                    "salmon_threads": 1,
                    "salmon_library_type": "A",
                    "salmon_extra_args": ["--validateMappings"],
                    "min_tpm_for_detected": 1.0,
                    "min_numreads_for_detected": 10.0,
                }
            }
        ),
        encoding="utf-8",
    )
    return config


def test_raw_read_intake_filters_to_retina_rnaseq_and_records_missing_pairs(tmp_path: Path):
    config_path = write_phase7d_fixture(tmp_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    phase = config["phase7d_rnaseq_quantification"]

    rows = build_raw_read_intake_rows(
        config,
        read_tsv(Path(phase["phase7a_manifest"]), PHASE7_RNASEQ_MANIFEST_COLUMNS),
        read_tsv(Path(phase["phase7a_readiness"]), PHASE7_RNASEQ_READINESS_COLUMNS),
    )

    assert [row["run"] for row in rows] == ["SRR1"]
    assert rows[0]["raw_pair_status"] == "PAIRED_FASTQ_MISSING_OR_INCOMPLETE"
    assert rows[0]["intake_status"] == "RAW_READS_NOT_AVAILABLE_LOCALLY"


def test_fastq_qc_reports_deterministic_lightweight_metrics(tmp_path: Path):
    config_path = write_phase7d_fixture(tmp_path, with_fastq=True)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    phase = config["phase7d_rnaseq_quantification"]
    intake = build_raw_read_intake_rows(
        config,
        read_tsv(Path(phase["phase7a_manifest"]), PHASE7_RNASEQ_MANIFEST_COLUMNS),
        read_tsv(Path(phase["phase7a_readiness"]), PHASE7_RNASEQ_READINESS_COLUMNS),
    )

    qc_rows = build_fastq_qc_rows(intake, max_records=10)
    by_read = {row["read"]: row for row in qc_rows}

    assert by_read["R1"]["records_inspected"] == "2"
    assert by_read["R1"]["mean_read_length"] == "4.00"
    assert by_read["R1"]["gc_percent"] == "50.000"
    assert by_read["R1"]["n_percent"] == "25.000"
    assert by_read["R1"]["qc_status"] == "PASS_LIGHTWEIGHT_FASTQ_QC"


def test_candidate_expression_matrix_uses_quant_sf_without_differential_claims(tmp_path: Path):
    config_path = write_phase7d_fixture(tmp_path, with_fastq=True)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    quant_dir = tmp_path / "salmon_quant" / "SRR1"
    quant_dir.mkdir(parents=True)
    quant_sf = quant_dir / "quant.sf"
    quant_sf.write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
        "tx1\t100\t80\t1.5\t10\n"
        "tx2\t100\t80\t2.5\t20\n"
        "tx3\t100\t80\t5.0\t50\n",
        encoding="utf-8",
    )

    rows = build_candidate_expression_matrix(
        config,
        [
            {
                "run": "SRR1",
                "tissue": "retina",
                "quant_sf": str(quant_sf),
                "quantification_status": "SALMON_QUANT_COMPLETE",
            }
        ],
    )
    by_gene = {row["gene_symbol"]: row for row in rows}

    assert by_gene["RAD51"]["summed_tpm"] == "4.000000"
    assert by_gene["RAD51"]["summed_numreads"] == "30.000"
    assert by_gene["RAD51"]["candidate_expression_status"] == "DETECTED_IN_RETINA_EXPLORATORY"
    assert by_gene["TP53"]["candidate_expression_status"] == "QUANTIFIED_REQUIRES_LOCUS_REVIEW_BEFORE_INTERPRETATION"
    assert "differential expression" in by_gene["RAD51"]["claim_language_guardrail"]


def test_salmon_quantification_reuses_existing_quant_sf_as_restart_checkpoint(tmp_path: Path):
    config_path = write_phase7d_fixture(tmp_path, with_fastq=True)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    phase = config["phase7d_rnaseq_quantification"]
    (Path(phase["salmon_index_dir"]) / "versionInfo.json").parent.mkdir(parents=True)
    (Path(phase["salmon_index_dir"]) / "versionInfo.json").write_text("{}", encoding="utf-8")
    quant_dir = Path(phase["salmon_quant_dir"]) / "SRR1"
    quant_dir.mkdir(parents=True)
    (quant_dir / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n", encoding="utf-8")
    intake = build_raw_read_intake_rows(
        config,
        read_tsv(Path(phase["phase7a_manifest"]), PHASE7_RNASEQ_MANIFEST_COLUMNS),
        read_tsv(Path(phase["phase7a_readiness"]), PHASE7_RNASEQ_READINESS_COLUMNS),
    )

    rows = run_salmon_quantification(config, intake, preflight_rows=[])

    assert rows[0]["quantification_status"] == "SALMON_QUANT_COMPLETE"
    assert rows[0]["salmon_command"] == "REUSED_EXISTING_QUANT_SF"
    assert "restart checkpoint" in rows[0]["message"]


def test_phase7d_cli_outputs_not_run_tables_when_reads_are_absent(tmp_path: Path):
    config_path = write_phase7d_fixture(tmp_path)
    outputs = {
        "intake": tmp_path / "intake.tsv",
        "fastq_qc": tmp_path / "fastq_qc.tsv",
        "preflight": tmp_path / "preflight.tsv",
        "quant": tmp_path / "quant.tsv",
        "matrix": tmp_path / "matrix.tsv",
        "report": tmp_path / "report.md",
    }

    run_phase7d(
        config_path,
        outputs["intake"],
        outputs["fastq_qc"],
        outputs["preflight"],
        outputs["quant"],
        outputs["matrix"],
        outputs["report"],
    )

    intake = read_tsv(outputs["intake"], PHASE7D_RAW_READ_INTAKE_COLUMNS)
    fastq_qc = read_tsv(outputs["fastq_qc"], PHASE7D_FASTQ_QC_COLUMNS)
    preflight = read_tsv(outputs["preflight"], PHASE7D_SALMON_PREFLIGHT_COLUMNS)
    quant = read_tsv(outputs["quant"], PHASE7D_SALMON_QUANT_SUMMARY_COLUMNS)
    matrix = read_tsv(outputs["matrix"], PHASE7D_CANDIDATE_EXPRESSION_MATRIX_COLUMNS)

    assert intake[0]["intake_status"] == "RAW_READS_NOT_AVAILABLE_LOCALLY"
    assert {row["qc_status"] for row in fastq_qc} == {"NOT_RUN_FASTQ_MISSING"}
    assert {row["check_id"]: row["status"] for row in preflight}["paired_fastqs"] == "MISSING"
    assert quant[0]["quantification_status"] == "NOT_RUN_MISSING_RAW_FASTQ_PAIR"
    assert {row["candidate_expression_status"] for row in matrix} == {"NOT_QUANTIFIED"}
    assert "retina-only" in outputs["report"].read_text(encoding="utf-8")
