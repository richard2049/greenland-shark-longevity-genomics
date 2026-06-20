from greenland_shark_longevity.local_file_inventory import detect_file_role, inventory_local_files


def test_detect_figshare_annotation_roles():
    role, detected_format, _notes = detect_file_role(type("P", (), {"name": "complete.proteins.faa", "suffix": ".faa"})())
    assert role == "genome_annotation_protein_fasta"
    assert detected_format == "protein_fasta"


def test_inventory_local_files_records_checksums_and_counts(tmp_path):
    protein = tmp_path / "complete.proteins.faa"
    protein.write_text(">p1\nMA\n>p2\nMKT\n", encoding="utf-8")
    rows = inventory_local_files(tmp_path, "RESOURCE", "SOURCE")
    assert len(rows) == 1
    assert rows[0]["file_name"] == "complete.proteins.faa"
    assert rows[0]["sequence_count"] == "2"
    assert len(rows[0]["md5"]) == 32
    assert len(rows[0]["sha256"]) == 64
