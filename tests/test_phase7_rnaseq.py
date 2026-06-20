from pathlib import Path

import yaml

from greenland_shark_longevity.phase7_rnaseq import (
    build_candidate_expression_plan,
    build_manifest_rows,
    build_readiness_rows,
    write_phase7_outputs,
)
from greenland_shark_longevity.schemas import (
    PHASE7_CANDIDATE_EXPRESSION_PLAN_COLUMNS,
    PHASE7_RNASEQ_MANIFEST_COLUMNS,
    PHASE7_RNASEQ_READINESS_COLUMNS,
)
from greenland_shark_longevity.utils import read_tsv


def write_panel(path: Path) -> None:
    path.write_text(
        yaml.safe_dump(
            {
                "panels": {
                    "dna_repair": {
                        "mechanism": "DNA repair/genome stability",
                        "candidates": [
                            {
                                "gene_symbol": "RAD51",
                                "synonyms": [],
                                "reference_ids": {"human_gene": "TODO"},
                                "caveats": "Expression support alone cannot establish repair activity.",
                            }
                        ],
                    },
                    "retina": {
                        "mechanism": "Tissue-specific preservation",
                        "candidates": [
                            {
                                "gene_symbol": "RH1",
                                "synonyms": ["rhodopsin"],
                                "reference_ids": {"human_gene": "TODO"},
                                "caveats": "Retinal expression remains tissue-specific.",
                            }
                        ],
                    },
                }
            }
        ),
        encoding="utf-8",
    )


def phase7_config(panel: Path) -> dict:
    return {
        "phase7_rnaseq": {
            "resource_id": "SMIC_RETINA_PRJNA1246101_2026",
            "bioproject": "PRJNA1246101",
            "sra_study": "SRP576249",
            "source_url": "https://www.ncbi.nlm.nih.gov/sra?from_uid=1246101&linkname=bioproject_sra_all",
            "retrieval_date": "2026-06-12",
            "raw_read_local_dir": "data/raw/rnaseq/SMIC_RETINA_PRJNA1246101_2026",
            "candidate_panel": str(panel),
            "minimum_replicates_for_expression_audit": 3,
            "registered_runs": [
                {
                    "run": "SRR1",
                    "experiment": "SRX1",
                    "sample": "SRS1",
                    "biosample": "SAMN1",
                    "sample_name": "Greenland_shark_retina_1",
                    "library_strategy": "RNA-Seq",
                    "library_source": "TRANSCRIPTOMIC",
                    "library_selection": "PolyA",
                    "library_layout": "PAIRED",
                    "platform": "ILLUMINA",
                    "model": "Illumina NovaSeq 6000",
                    "spots": 10,
                    "bases": 2000,
                    "size_mb": 1,
                    "taxon_id": "191813",
                    "scientific_name": "Somniosus microcephalus",
                    "sex": "female",
                    "tissue": "retina",
                },
                {
                    "run": "SRR2",
                    "experiment": "SRX2",
                    "sample": "SRS2",
                    "biosample": "SAMN2",
                    "sample_name": "Greenland_shark_retina_2",
                    "library_strategy": "RNA-Seq",
                    "library_source": "TRANSCRIPTOMIC",
                    "library_selection": "PolyA",
                    "library_layout": "PAIRED",
                    "platform": "ILLUMINA",
                    "model": "Illumina NovaSeq 6000",
                    "spots": 10,
                    "bases": 2000,
                    "size_mb": 1,
                    "taxon_id": "191813",
                    "scientific_name": "Somniosus microcephalus",
                    "sex": "male",
                    "tissue": "retina",
                },
                {
                    "run": "SRR3",
                    "experiment": "SRX3",
                    "sample": "SRS3",
                    "biosample": "SAMN3",
                    "sample_name": "Greenland_shark_retina_3",
                    "library_strategy": "RNA-Seq",
                    "library_source": "TRANSCRIPTOMIC",
                    "library_selection": "PolyA",
                    "library_layout": "PAIRED",
                    "platform": "ILLUMINA",
                    "model": "Illumina NovaSeq 6000",
                    "spots": 10,
                    "bases": 2000,
                    "size_mb": 1,
                    "taxon_id": "191813",
                    "scientific_name": "Somniosus microcephalus",
                    "sex": "male",
                    "tissue": "retina",
                },
                {
                    "run": "SRR4",
                    "experiment": "SRX4",
                    "sample": "SRS4",
                    "biosample": "SAMN4",
                    "sample_name": "Greenland_shark_retina_wgs",
                    "library_strategy": "WGS",
                    "library_source": "GENOMIC",
                    "library_selection": "other",
                    "library_layout": "PAIRED",
                    "platform": "ILLUMINA",
                    "model": "Illumina NovaSeq 6000",
                    "spots": 10,
                    "bases": 2000,
                    "size_mb": 1,
                    "taxon_id": "191813",
                    "scientific_name": "Somniosus microcephalus",
                    "sex": "unknown",
                    "tissue": "retina",
                },
            ],
        }
    }


def test_phase7_runinfo_metadata_classification(tmp_path: Path):
    panel = tmp_path / "candidate_panels.yaml"
    write_panel(panel)
    config = phase7_config(panel)

    manifest_rows = build_manifest_rows(config)
    readiness_rows = build_readiness_rows(manifest_rows)
    by_run = {row["run"]: row for row in readiness_rows}

    assert len(manifest_rows) == 4
    assert by_run["SRR1"]["expression_readiness_status"] == "RNA_SEQ_METADATA_READY_FOR_CANDIDATE_EXPRESSION_AUDIT"
    assert by_run["SRR4"]["expression_readiness_status"] == "EXCLUDED_NON_TRANSCRIPTOMIC_RESOURCE"
    assert by_run["SRR4"]["included_for_candidate_expression_plan"] == "false"
    assert "quantify expression" in by_run["SRR1"]["conservative_interpretation"]


def test_phase7_candidate_plan_is_metadata_only(tmp_path: Path):
    panel = tmp_path / "candidate_panels.yaml"
    write_panel(panel)
    config = phase7_config(panel)
    readiness_rows = build_readiness_rows(build_manifest_rows(config))

    plan_rows = build_candidate_expression_plan(panel, readiness_rows, min_replicates=3)
    by_gene = {row["gene_symbol"]: row for row in plan_rows}

    assert by_gene["RAD51"]["usable_rnaseq_run_count"] == "3"
    assert by_gene["RAD51"]["usable_tissues"] == "retina"
    assert by_gene["RH1"]["current_expression_status"] == "NOT_QUANTIFIED_PHASE7A_METADATA_ONLY"
    assert "activated" not in by_gene["RH1"]["conservative_interpretation"].lower()


def test_phase7_cli_outputs_expected_tables(tmp_path: Path):
    panel = tmp_path / "candidate_panels.yaml"
    config_path = tmp_path / "config.yaml"
    manifest = tmp_path / "rnaseq_manifest.tsv"
    readiness = tmp_path / "readiness.tsv"
    plan = tmp_path / "plan.tsv"
    report = tmp_path / "report.md"
    write_panel(panel)
    config_path.write_text(yaml.safe_dump(phase7_config(panel)), encoding="utf-8")

    write_phase7_outputs(config_path, manifest, readiness, plan, report)

    assert len(read_tsv(manifest, PHASE7_RNASEQ_MANIFEST_COLUMNS)) == 4
    assert len(read_tsv(readiness, PHASE7_RNASEQ_READINESS_COLUMNS)) == 4
    assert len(read_tsv(plan, PHASE7_CANDIDATE_EXPRESSION_PLAN_COLUMNS)) == 2
    assert "does not download raw reads" in report.read_text(encoding="utf-8")
