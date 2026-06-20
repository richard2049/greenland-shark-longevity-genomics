from pathlib import Path

import yaml

from greenland_shark_longevity.phase6_telomere import (
    build_motif_scan_rows,
    build_telomere_gene_audit,
    summarize_enrichment,
    write_phase6_outputs,
)
from greenland_shark_longevity.schemas import (
    COPY_NUMBER_COLUMNS,
    DUPLICATION_AUDIT_COLUMNS,
    EVIDENCE_COLUMNS,
    PHASE6_SCAFFOLD_END_ENRICHMENT_COLUMNS,
    PHASE6_TELOMERE_GENE_AUDIT_COLUMNS,
    PHASE6_TELOMERE_MOTIF_SCAN_COLUMNS,
)
from greenland_shark_longevity.utils import read_tsv, write_tsv


def write_panel(path: Path) -> None:
    path.write_text(
        yaml.safe_dump(
            {
                "panels": {
                    "telomere_shelterin": {
                        "mechanism": "Telomere-related biology",
                        "candidates": [
                            {
                                "gene_symbol": "TERT",
                                "synonyms": ["telomerase reverse transcriptase"],
                                "reference_ids": {"human_gene": "TODO", "uniprot": "TODO"},
                                "caveats": "Presence does not imply telomerase activity or telomere length.",
                            },
                            {
                                "gene_symbol": "TERF2",
                                "synonyms": ["TRF2"],
                                "reference_ids": {"human_gene": "TODO", "uniprot": "TODO"},
                                "caveats": "Panel-only fixture.",
                            },
                        ],
                    }
                }
            }
        ),
        encoding="utf-8",
    )


def test_telomeric_motif_scan_counts_terminal_and_control_windows(tmp_path: Path):
    fasta = tmp_path / "genome.fna"
    fasta.write_text(
        ">scaffold_1\n"
        "TTAGGGTTAGGG" + "A" * 20 + "CCCTAA" + "C" * 20 + "CCCTAACCCTAA\n"
        ">short\n"
        "TTAGGG\n",
        encoding="utf-8",
    )

    rows = build_motif_scan_rows(fasta, "TEST", ["TTAGGG", "CCCTAA"], terminal_window_bp=12, min_sequence_length_bp=10)
    summary = summarize_enrichment(rows, "TEST", fasta, ["TTAGGG", "CCCTAA"], 12)

    assert len(rows) == 1
    assert rows[0]["left_TTAGGG_count"] == "2"
    assert rows[0]["right_CCCTAA_count"] == "2"
    assert rows[0]["left_control_CCCTAA_count"] == "0"
    assert summary[0]["terminal_total_motif_count"] == "4"
    assert summary[0]["motif_scan_status"] == "CANONICAL_MOTIF_END_ENRICHMENT_RECORDED_SEQUENCE_CONTEXT_ONLY"


def test_telomere_gene_audit_distinguishes_assessed_and_panel_only_genes(tmp_path: Path):
    panel = tmp_path / "candidate_panels.yaml"
    copy = tmp_path / "copy.tsv"
    audit = tmp_path / "audit.tsv"
    evidence = tmp_path / "evidence.tsv"
    write_panel(panel)
    write_tsv(
        copy,
        [
            {
                "mechanism": "Telomere-related biology",
                "gene_symbol": "TERT",
                "species_id": "smic",
                "resource_id": "SMIC",
                "orthogroup_id": "OG1",
                "copy_count": "1",
                "mapping_status": "ANNOTATION_SYMBOL_MATCH_ORTHOGROUP_MAPPED",
                "resources_supporting": "SMIC",
            }
        ],
        COPY_NUMBER_COLUMNS,
    )
    write_tsv(
        audit,
        [
            {
                "mechanism": "Telomere-related biology",
                "gene_symbol": "TERT",
                "orthogroup_id": "OG1",
                "copy_count": "1",
                "artifact_risk": "low",
            }
        ],
        DUPLICATION_AUDIT_COLUMNS,
    )
    write_tsv(
        evidence,
        [
            {
                "mechanism": "Telomere-related biology",
                "gene_or_pathway": "TERT",
                "evidence_type": "candidate_presence_or_not_assessed",
                "evidence_tier": "Tier 2",
                "resources_supporting": "SMIC",
                "artifact_risk": "low",
                "biological_interpretation": "Plausible but incomplete computational evidence.",
                "relevance_to_aging_longevity": "Hypothesis-generating; requires validation before biological interpretation.",
                "translational_relevance": "NOT_ASSESSED",
                "required_validation": "REQUIRES_VALIDATION",
                "claim_language_guardrail": "test",
            }
        ],
        EVIDENCE_COLUMNS,
    )

    rows = build_telomere_gene_audit(panel, copy, audit, evidence)
    by_gene = {row["gene_symbol"]: row for row in rows}

    assert by_gene["TERT"]["readiness_status"] == "ORTHOLOGY_EVIDENCE_READY_FOR_DOMAIN_AND_CONTEXT_VALIDATION"
    assert by_gene["TERF2"]["readiness_status"] == "PANEL_ONLY_NOT_ASSESSED_IN_CURRENT_ORTHOLOGY_OUTPUT"
    assert "absence" not in by_gene["TERF2"]["conservative_interpretation"].lower()


def test_phase6_cli_outputs_expected_tables(tmp_path: Path):
    panel = tmp_path / "candidate_panels.yaml"
    copy = tmp_path / "copy.tsv"
    audit = tmp_path / "audit.tsv"
    evidence = tmp_path / "evidence.tsv"
    fasta = tmp_path / "genome.fna"
    config = tmp_path / "config.yaml"
    write_panel(panel)
    fasta.write_text(">scaffold_1\nTTAGGGTTAGGG" + "A" * 60 + "CCCTAACCCTAA\n", encoding="utf-8")
    write_tsv(copy, [], COPY_NUMBER_COLUMNS)
    write_tsv(audit, [], DUPLICATION_AUDIT_COLUMNS)
    write_tsv(evidence, [], EVIDENCE_COLUMNS)
    motif_out = tmp_path / "motif.tsv"
    enrich_out = tmp_path / "enrich.tsv"
    gene_out = tmp_path / "genes.tsv"
    report = tmp_path / "report.md"
    config.write_text(
        yaml.safe_dump(
            {
                "phase6_telomere": {
                    "genome_fasta": str(fasta),
                    "resource_id": "TEST",
                    "motifs": ["TTAGGG", "CCCTAA"],
                    "terminal_window_bp": 12,
                    "min_sequence_length_bp": 10,
                    "candidate_panel": str(panel),
                    "candidate_copy_number": str(copy),
                    "candidate_duplication_audit": str(audit),
                    "integrated_evidence": str(evidence),
                }
            }
        ),
        encoding="utf-8",
    )

    write_phase6_outputs(config, motif_out, enrich_out, gene_out, report)

    assert read_tsv(motif_out, PHASE6_TELOMERE_MOTIF_SCAN_COLUMNS)
    assert read_tsv(enrich_out, PHASE6_SCAFFOLD_END_ENRICHMENT_COLUMNS)
    assert len(read_tsv(gene_out, PHASE6_TELOMERE_GENE_AUDIT_COLUMNS)) == 2
    assert "does not infer telomere length" in report.read_text(encoding="utf-8")
