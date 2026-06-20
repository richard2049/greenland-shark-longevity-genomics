from pathlib import Path

from greenland_shark_longevity.fasta import assembly_qc, n50, parse_fasta, protein_qc


def test_n50_calculation():
    assert n50([50, 28, 12]) == 50
    assert n50([10, 10, 8]) == 10


def test_assembly_qc_demo_fixture():
    row = assembly_qc(Path("data/demo/assemblies/demo_smic_assembly.fa"), "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE")
    assert row["sequence_count"] == 3
    assert row["total_length_bp"] == 90
    assert row["n50_bp"] == 50
    assert row["n_count"] == 6
    assert row["busco_complete_percent"] == "NOT_ASSESSED"


def test_duplicate_fasta_ids_are_reported(tmp_path):
    fasta = tmp_path / "duplicate.fa"
    fasta.write_text(">a\nACG\n>a\nTTT\n", encoding="utf-8")
    records, duplicates = parse_fasta(fasta)
    assert len(records) == 2
    assert duplicates == ["a"]


def test_protein_qc_reports_ambiguous_and_stop():
    row = protein_qc(Path("data/demo/proteins/demo_smic_proteins.fa"), "DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE")
    assert row["protein_count"] == 4
    assert row["ambiguous_residue_count"] == 2
    assert row["sequences_with_stop_codon"] == 1

