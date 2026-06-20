from pathlib import Path

from greenland_shark_longevity.phase7e_expression_hardening import (
    build_candidate_hardened_rows,
    build_parameter_review_rows,
    build_run_qc_rows,
)


def phase7e_config(tmp_path: Path) -> dict:
    return {
        "phase7d_rnaseq_quantification": {
            "salmon_executable": "tools/salmon-docker.cmd",
            "salmon_index_dir": str(tmp_path / "salmon_index"),
            "salmon_extra_args": ["--validateMappings", "--seqBias", "--gcBias"],
            "salmon_library_type": "A",
            "min_tpm_for_detected": 1.0,
            "min_numreads_for_detected": 10.0,
            "fastq_qc_max_records": 100000,
        },
        "phase7e_expression_hardening": {
            "low_mapping_rate_warning_percent": 30.0,
            "very_low_mapping_rate_percent": 20.0,
            "min_detected_runs_for_consistent_retina_support": 2,
            "min_detection_fraction_for_consistent_retina_support": 0.66,
        },
    }


def write_salmon_meta(tmp_path: Path, run: str, percent_mapped: float, num_processed: int = 1000) -> dict[str, str]:
    output_dir = tmp_path / "salmon" / run
    aux = output_dir / "aux_info"
    aux.mkdir(parents=True)
    (aux / "meta_info.json").write_text(
        (
            "{"
            f'"num_processed": {num_processed}, '
            f'"num_mapped": {int(num_processed * percent_mapped / 100)}, '
            f'"percent_mapped": {percent_mapped}, '
            '"library_types": ["ISR"]'
            "}"
        ),
        encoding="utf-8",
    )
    stderr = tmp_path / f"{run}.stderr.log"
    stderr.write_text("[jointLog] [info] Number of decoys : 0\n", encoding="utf-8")
    return {
        "run": run,
        "tissue": "retina",
        "read1_path": f"{run}_1.fastq.gz",
        "read2_path": f"{run}_2.fastq.gz",
        "salmon_output_dir": str(output_dir),
        "quant_sf": str(output_dir / "quant.sf"),
        "quantification_status": "SALMON_QUANT_COMPLETE",
        "salmon_command": "salmon quant",
        "log_stdout": str(tmp_path / f"{run}.stdout.log"),
        "log_stderr": str(stderr),
        "message": "complete",
    }


def fastq_qc_rows(run: str) -> list[dict[str, str]]:
    rows = []
    for read in ["R1", "R2"]:
        rows.append(
            {
                "run": run,
                "read": read,
                "fastq_path": f"{run}_{read}.fastq.gz",
                "file_status": "LOCAL_AVAILABLE",
                "records_inspected": "100000",
                "mean_read_length": "100.0",
                "gc_percent": "45.0",
                "n_percent": "0.001",
                "mean_phred_quality": "36.0",
                "qc_status": "PASS_LIGHTWEIGHT_FASTQ_QC",
                "notes": "fixture",
            }
        )
    return rows


def test_run_qc_flags_mapping_rates_and_no_decoy_index(tmp_path: Path):
    config = phase7e_config(tmp_path)
    quant_rows = [
        write_salmon_meta(tmp_path, "SRR_PASS", 32.8),
        write_salmon_meta(tmp_path, "SRR_LOW", 27.8),
        write_salmon_meta(tmp_path, "SRR_VERY_LOW", 14.4),
    ]
    qc_rows = fastq_qc_rows("SRR_PASS") + fastq_qc_rows("SRR_LOW") + fastq_qc_rows("SRR_VERY_LOW")

    rows = build_run_qc_rows(config, quant_rows, qc_rows)
    by_run = {row["run"]: row for row in rows}

    assert by_run["SRR_PASS"]["mapping_rate_status"] == "PASS_FIRST_PASS_MAPPING_RATE"
    assert by_run["SRR_LOW"]["mapping_rate_status"] == "LOW_MAPPING_RATE_CAUTIOUS_USE"
    assert by_run["SRR_VERY_LOW"]["mapping_rate_status"] == "VERY_LOW_MAPPING_RATE_REVIEW_REQUIRED"
    assert {row["salmon_decoy_count"] for row in rows} == {"0"}

    parameter_rows = build_parameter_review_rows(config, rows)
    decoy_row = {row["parameter"]: row for row in parameter_rows}["salmon_index_decoys"]
    assert decoy_row["review_status"] == "FIRST_PASS_ACCEPTABLE_CONFIRMATORY_DECOY_INDEX_RECOMMENDED"


def test_candidate_hardening_allows_low_ambiguity_support_but_defers_locus_ambiguous_gene(tmp_path: Path):
    config = phase7e_config(tmp_path)
    run_qc_rows = [
        {
            "run": "SRR1",
            "salmon_percent_mapped": "32.0",
            "mapping_rate_status": "PASS_FIRST_PASS_MAPPING_RATE",
        },
        {
            "run": "SRR2",
            "salmon_percent_mapped": "28.0",
            "mapping_rate_status": "LOW_MAPPING_RATE_CAUTIOUS_USE",
        },
        {
            "run": "SRR3",
            "salmon_percent_mapped": "14.0",
            "mapping_rate_status": "VERY_LOW_MAPPING_RATE_REVIEW_REQUIRED",
        },
    ]
    expression_rows = [
        {
            "mechanism": "DNA repair/genome stability",
            "gene_symbol": "ERCC4",
            "run": "SRR1",
            "tissue": "retina",
            "phase7c_quantification_readiness": "REFERENCE_READY_FOR_FIRST_PASS_QUANTIFICATION",
            "matched_transcript_count": "1",
            "quantified_transcript_count": "1",
            "summed_tpm": "3.0",
            "summed_numreads": "100",
            "candidate_expression_status": "DETECTED_IN_RETINA_EXPLORATORY",
            "quantification_status": "SALMON_QUANT_COMPLETE",
            "supporting_files": "quant.sf",
        },
        {
            "mechanism": "DNA repair/genome stability",
            "gene_symbol": "ERCC4",
            "run": "SRR2",
            "tissue": "retina",
            "phase7c_quantification_readiness": "REFERENCE_READY_FOR_FIRST_PASS_QUANTIFICATION",
            "matched_transcript_count": "1",
            "quantified_transcript_count": "1",
            "summed_tpm": "2.0",
            "summed_numreads": "80",
            "candidate_expression_status": "DETECTED_IN_RETINA_EXPLORATORY",
            "quantification_status": "SALMON_QUANT_COMPLETE",
            "supporting_files": "quant.sf",
        },
        {
            "mechanism": "DNA repair/genome stability",
            "gene_symbol": "ERCC4",
            "run": "SRR3",
            "tissue": "retina",
            "phase7c_quantification_readiness": "REFERENCE_READY_FOR_FIRST_PASS_QUANTIFICATION",
            "matched_transcript_count": "1",
            "quantified_transcript_count": "1",
            "summed_tpm": "0.0",
            "summed_numreads": "0",
            "candidate_expression_status": "NOT_DETECTED_UNDER_CURRENT_RETINA_QUANTIFICATION",
            "quantification_status": "SALMON_QUANT_COMPLETE",
            "supporting_files": "quant.sf",
        },
        {
            "mechanism": "DNA repair/genome stability",
            "gene_symbol": "RAD51",
            "run": "SRR1",
            "tissue": "retina",
            "phase7c_quantification_readiness": "REFERENCE_PRESENT_TARGETED_LOCUS_VALIDATION_REQUIRED",
            "matched_transcript_count": "4",
            "quantified_transcript_count": "4",
            "summed_tpm": "10.0",
            "summed_numreads": "120",
            "candidate_expression_status": "QUANTIFIED_REQUIRES_LOCUS_REVIEW_BEFORE_INTERPRETATION",
            "quantification_status": "SALMON_QUANT_COMPLETE",
            "supporting_files": "quant.sf",
        },
    ]
    phase7c_rows = [
        {
            "gene_symbol": "ERCC4",
            "phase7c_quantification_readiness": "REFERENCE_READY_FOR_FIRST_PASS_QUANTIFICATION",
            "reference_mapping_status": "TRANSCRIPT_REFERENCE_PRESENT",
            "ambiguity_status": "low_reference_mapping_risk",
            "matched_transcript_count": "1",
            "required_validation": "review",
        },
        {
            "gene_symbol": "RAD51",
            "phase7c_quantification_readiness": "REFERENCE_PRESENT_TARGETED_LOCUS_VALIDATION_REQUIRED",
            "reference_mapping_status": "TRANSCRIPT_REFERENCE_PRESENT_TARGETED_REVIEW_REQUIRED",
            "ambiguity_status": "high_locus_or_paralog_ambiguity",
            "matched_transcript_count": "4",
            "required_validation": "targeted review",
        },
    ]
    duplication_rows = [
        {
            "gene_symbol": "ERCC4",
            "mapping_status": "ANNOTATION_SYMBOL_MATCH_ORTHOGROUP_MAPPED",
            "artifact_risk": "low",
            "required_validation": "orthology review",
        },
        {
            "gene_symbol": "RAD51",
            "mapping_status": "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH",
            "artifact_risk": "not_assessable",
            "required_validation": "manual review",
        },
    ]
    phase4e_rows = [
        {
            "gene_symbol": "RAD51",
            "phase4e_hardened_status": "HARDENED_CANDIDATE_LOCI_REQUIRE_CROSS_RESOURCE_VALIDATION",
            "artifact_risk": "moderate",
            "required_validation": "locus review",
            "supporting_files": "phase4e.tsv",
        }
    ]

    rows = build_candidate_hardened_rows(config, expression_rows, phase7c_rows, duplication_rows, phase4e_rows, [], run_qc_rows)
    by_gene = {row["gene_symbol"]: row for row in rows}

    assert by_gene["ERCC4"]["expression_interpretation_status"] == "RETINA_DETECTED_FIRST_PASS_WITH_QC_CAVEATS"
    assert by_gene["ERCC4"]["phase7e_expression_support_level"] == "CAUTIOUS_RETINA_EXPRESSION_SUPPORT_FOR_PHASE8"
    assert "very_low_mapping_rate_caveat" in by_gene["ERCC4"]["artifact_risk"]
    assert by_gene["RAD51"]["expression_interpretation_status"] == "LOCUS_REVIEW_REQUIRED_BEFORE_EXPRESSION_SUPPORT"
    assert by_gene["RAD51"]["phase7e_expression_support_level"] == "DEFER_PHASE8_EXPRESSION_SUPPORT_UNTIL_LOCUS_REVIEW"
