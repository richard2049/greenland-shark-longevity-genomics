from pathlib import Path

from greenland_shark_longevity.phase5b_repeat_annotation import (
    build_output_inventory_rows,
    build_plan_rows,
    build_preflight_rows,
)


def minimal_config(tmp_path: Path, genome_exists: bool = True) -> dict:
    genome = tmp_path / "genome.fna.gz"
    if genome_exists:
        genome.write_bytes(b"not really gzip but present")
    return {
        "phase5b_repeat_annotation": {
            "genome_fasta": str(genome),
            "work_dir": str(tmp_path / "work"),
            "output_dir": str(tmp_path / "out"),
            "repeatmasker_output_dir": str(tmp_path / "out" / "repeatmasker_raw"),
            "staged_genome_fasta": str(tmp_path / "work" / "genome.fna"),
            "database_name": "smic_test",
            "container_image": "dfam/tetools:latest",
            "threads": 2,
            "custom_library": str(tmp_path / "out" / "families.fa"),
            "repeatmasker_gff": str(tmp_path / "out" / "repeatmasker.out.gff"),
            "repeatmasker_out": str(tmp_path / "out" / "repeatmasker.out"),
            "repeatmasker_tbl": str(tmp_path / "out" / "repeatmasker.tbl"),
        }
    }


def test_phase5b_preflight_distinguishes_present_genome_and_missing_outputs(tmp_path: Path):
    rows = build_preflight_rows(minimal_config(tmp_path))
    by_check = {row["check_id"]: row for row in rows}

    assert by_check["PHASE5B-GENOME-FASTA"]["status"] == "PASS"
    assert by_check["PHASE5B-REPEATMASKER-GFF"]["status"] == "EXPECTED_NOT_PRESENT"
    assert by_check["PHASE5B-CONTAINER-IMAGE"]["status"] == "WARN"


def test_phase5b_plan_is_manual_and_contains_repeatmodeler_repeatmasker_steps(tmp_path: Path):
    rows = build_plan_rows(minimal_config(tmp_path))
    text = "\n".join(row["command_template"] for row in rows)

    assert "BuildDatabase" in text
    assert "RepeatModeler" in text
    assert "RepeatMasker" in text
    assert all(row["run_status"] != "RUN_AUTOMATICALLY" for row in rows)


def test_phase5b_inventory_marks_repeatmasker_gff_as_phase5_import_candidate(tmp_path: Path):
    config = minimal_config(tmp_path)
    gff = Path(config["phase5b_repeat_annotation"]["repeatmasker_gff"])
    gff.parent.mkdir(parents=True)
    gff.write_text("##gff-version 3\n", encoding="utf-8")

    rows = build_output_inventory_rows(config)
    by_id = {row["output_id"]: row for row in rows}

    assert by_id["repeatmasker_gff"]["status"] == "LOCAL_AVAILABLE"
    assert by_id["repeatmasker_gff"]["recommended_for_phase5_import"] == "true"
    assert len(by_id["repeatmasker_gff"]["sha256"]) == 64
