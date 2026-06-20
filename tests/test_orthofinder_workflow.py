from pathlib import Path

from greenland_shark_longevity.orthofinder_workflow import (
    build_preflight_rows,
    missing_orthofinder_message,
    parse_gene_count_matrix,
)
from greenland_shark_longevity.schemas import ORTHOFINDER_INPUT_MANIFEST_COLUMNS
from greenland_shark_longevity.utils import write_tsv


def test_preflight_counts_ready_species(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "orthofinder:",
                "  executable: definitely_missing_orthofinder_binary",
                "  minimum_species: 2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    fasta_a = tmp_path / "a.faa"
    fasta_b = tmp_path / "b.faa"
    fasta_a.write_text(">a\nMA\n", encoding="utf-8")
    fasta_b.write_text(">b\nMA\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    write_tsv(
        manifest,
        [
            {
                "species_id": "a",
                "resource_id": "A",
                "source_protein_fasta": "src_a",
                "orthofinder_input_path": str(fasta_a),
                "protein_set_scope": "GENOME_WIDE_PROTEIN_FASTA_CANDIDATE",
                "protein_count": "1",
                "qc_status": "LOCAL_AVAILABLE",
                "orthofinder_ready": "True",
                "blocker": "NONE",
                "next_action": "run",
                "notes": "",
            },
            {
                "species_id": "b",
                "resource_id": "B",
                "source_protein_fasta": "src_b",
                "orthofinder_input_path": str(fasta_b),
                "protein_set_scope": "GENOME_WIDE_PROTEIN_FASTA_CANDIDATE",
                "protein_count": "1",
                "qc_status": "LOCAL_AVAILABLE",
                "orthofinder_ready": "True",
                "blocker": "NONE",
                "next_action": "run",
                "notes": "",
            },
        ],
        ORTHOFINDER_INPUT_MANIFEST_COLUMNS,
    )
    rows = build_preflight_rows(config, manifest)
    by_check = {row["check_name"]: row for row in rows}
    assert by_check["orthofinder_binary"]["status"] == "FAIL"
    assert by_check["ready_species_count"]["status"] == "PASS"


def test_missing_orthofinder_message_explains_native_windows_blocker(monkeypatch):
    monkeypatch.setattr("greenland_shark_longevity.orthofinder_workflow.platform.system", lambda: "Windows")
    message = missing_orthofinder_message("orthofinder")
    assert "native Windows" in message
    assert "WSL/Linux" in message


def test_parse_gene_count_matrix(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    write_tsv(
        manifest,
        [
            {
                "species_id": "smic",
                "resource_id": "SMIC",
                "source_protein_fasta": "src_smic",
                "orthofinder_input_path": str(tmp_path / "smic__SMIC.faa"),
                "protein_set_scope": "GENOME_WIDE_PROTEIN_FASTA_CANDIDATE",
                "protein_count": "2",
                "qc_status": "LOCAL_AVAILABLE",
                "orthofinder_ready": "True",
                "blocker": "NONE",
                "next_action": "run",
                "notes": "",
            },
            {
                "species_id": "cmil",
                "resource_id": "CMIL",
                "source_protein_fasta": "src_cmil",
                "orthofinder_input_path": str(tmp_path / "cmil__CMIL.faa"),
                "protein_set_scope": "GENOME_WIDE_PROTEIN_FASTA_CANDIDATE",
                "protein_count": "2",
                "qc_status": "LOCAL_AVAILABLE",
                "orthofinder_ready": "True",
                "blocker": "NONE",
                "next_action": "run",
                "notes": "",
            },
        ],
        ORTHOFINDER_INPUT_MANIFEST_COLUMNS,
    )
    orthogroups = tmp_path / "Results_Test" / "Orthogroups"
    orthogroups.mkdir(parents=True)
    (orthogroups / "Orthogroups.GeneCount.tsv").write_text(
        "Orthogroup\tsmic__SMIC\tcmil__CMIL\tTotal\n"
        "OG0000001\t2\t1\t3\n"
        "OG0000002\t0\t2\t2\n",
        encoding="utf-8",
    )
    long_rows, summary_rows, result_root = parse_gene_count_matrix(tmp_path, manifest)
    assert result_root.name == "Results_Test"
    assert len(long_rows) == 4
    summary_by_species = {row["species_id"]: row for row in summary_rows}
    assert summary_by_species["smic"]["orthogroups_with_genes"] == "1"
    assert summary_by_species["cmil"]["total_gene_copies_in_orthogroups"] == "3"
