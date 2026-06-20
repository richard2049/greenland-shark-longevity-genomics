from pathlib import Path

import yaml

from greenland_shark_longevity.phase5c_repeat_qc import parse_repeatmasker_out_features, write_phase5c_outputs
from greenland_shark_longevity.phase5_repeat_context import CandidateWindow, load_assembly_seqid_aliases
from greenland_shark_longevity.schemas import (
    PHASE4E_LOCUS_HARDENING_COLUMNS,
    PHASE5B_OUTPUT_INVENTORY_COLUMNS,
    PHASE5C_GENE_QC_COLUMNS,
    PHASE5C_INTEGRITY_COLUMNS,
    PHASE5C_LOCUS_QC_COLUMNS,
    PHASE5_LOCUS_REPEAT_CONTEXT_COLUMNS,
)
from greenland_shark_longevity.utils import read_tsv, write_tsv


def write_repeatmasker_out(path: Path) -> None:
    path.write_text(
        "   SW   perc perc perc  query              position in query                 matching           repeat                position in repeat\n"
        " score   div. del. ins.  sequence           begin     end            (left)   repeat             class/family      begin   end    (left)      ID\n"
        "\n"
        " 1000    1.0  0.0  0.0  JBLTJD010000033.1        120      180      (0) + rnd-1_family-1     LINE/CR1        1 61 (0) 1\n"
        "  500    2.0  0.0  0.0  JBLTJD010000033.1        230      260      (0) + rnd-2_family-1     Unknown        1 31 (0) 2\n"
        "  900    1.0  0.0  0.0  JBLTJD010000012.1        500      600      (0) + rnd-3_family-1     DNA/hAT        1 101 (0) 3\n",
        encoding="utf-8",
    )


def write_candidate_loci(path: Path) -> None:
    write_tsv(
        path,
        [
            {
                "gene_symbol": "FTH1B",
                "mechanism": "Ferroptosis, iron handling, and oxidative stress",
                "locus_cluster_id": "FTH1B_L001",
                "target_seqid": "JBLTJD010000033.1",
                "annotation_seqid": "scaffold_33",
                "locus_start": "100",
                "locus_end": "200",
                "locus_span_bp": "101",
            }
        ],
        PHASE4E_LOCUS_HARDENING_COLUMNS,
    )


def write_phase5_locus_context(path: Path) -> None:
    write_tsv(
        path,
        [
            {
                "gene_symbol": "FTH1B",
                "mechanism": "Ferroptosis, iron handling, and oxidative stress",
                "locus_cluster_id": "FTH1B_L001",
                "annotation_seqid": "scaffold_33",
                "locus_start": "100",
                "locus_end": "200",
                "locus_span_bp": "101",
                "window_start": "90",
                "window_end": "210",
                "window_size_bp": "121",
                "repeat_annotation_status": "REPEAT_ANNOTATION_AVAILABLE",
                "repeat_overlap_count": "1",
                "repeat_overlap_bp": "61",
                "repeat_overlap_fraction": "0.60396",
            }
        ],
        PHASE5_LOCUS_REPEAT_CONTEXT_COLUMNS,
    )


def test_repeatmasker_out_parser_normalizes_accessions_and_filters_windows(tmp_path: Path):
    out = tmp_path / "smic_tokyo.repeatmasker.out"
    report = tmp_path / "assembly_report.txt"
    write_repeatmasker_out(out)
    report.write_text(
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
        "scaffold_33\tunplaced-scaffold\tna\tna\tJBLTJD010000033.1\t<>\tna\tPrimary Assembly\t1000\tna\n",
        encoding="utf-8",
    )
    aliases = load_assembly_seqid_aliases(report)
    windows = {"scaffold_33": [CandidateWindow("scaffold_33", 90, 210)]}

    features = parse_repeatmasker_out_features(out, aliases, windows)

    assert len(features) == 1
    assert features[0].seqid == "scaffold_33"
    assert features[0].repeat_class == "LINE"
    assert features[0].repeat_family == "CR1"


def test_phase5c_writes_integrity_locus_and_gene_qc(tmp_path: Path):
    candidate_loci = tmp_path / "phase4e_loci.tsv"
    phase5_locus = tmp_path / "phase5_locus.tsv"
    phase5b_inventory = tmp_path / "phase5b_inventory.tsv"
    out = tmp_path / "smic_tokyo.repeatmasker.out"
    gff = tmp_path / "smic_tokyo.repeatmasker.out.gff"
    families = tmp_path / "smic_tokyo-families.fa"
    log = tmp_path / "repeatmasker.log"
    config = tmp_path / "config.yaml"
    assembly_report = tmp_path / "assembly_report.txt"

    write_candidate_loci(candidate_loci)
    write_phase5_locus_context(phase5_locus)
    write_repeatmasker_out(out)
    gff.write_text("##gff-version 3\n", encoding="utf-8")
    families.write_text(">family1\nACGT\n", encoding="utf-8")
    log.write_text("forksys: Program terminated by a signal 9.\n", encoding="utf-8")
    assembly_report.write_text(
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
        "scaffold_33\tunplaced-scaffold\tna\tna\tJBLTJD010000033.1\t<>\tna\tPrimary Assembly\t1000\tna\n",
        encoding="utf-8",
    )
    write_tsv(
        phase5b_inventory,
        [
            {
                "output_id": "repeatmasker_out",
                "file_role": "repeatmasker_raw_out",
                "local_path": str(out),
                "status": "PRESENT",
                "byte_size": str(out.stat().st_size),
                "sha256": "TEST_SHA256",
            }
        ],
        PHASE5B_OUTPUT_INVENTORY_COLUMNS,
    )
    config.write_text(
        yaml.safe_dump(
            {
                "phase5_repeat_context": {
                    "assembly_report": str(assembly_report),
                    "candidate_window_bp": 10,
                },
                "phase5b_repeat_annotation": {
                    "custom_library": str(families),
                    "repeatmasker_gff": str(gff),
                    "repeatmasker_out": str(out),
                    "repeatmasker_tbl": str(tmp_path / "missing.tbl"),
                },
                "phase5c_repeat_qc": {
                    "log_files": [str(log)],
                },
            }
        ),
        encoding="utf-8",
    )

    integrity_output = tmp_path / "integrity.tsv"
    locus_output = tmp_path / "locus_qc.tsv"
    gene_output = tmp_path / "gene_qc.tsv"
    report_output = tmp_path / "report.md"
    write_phase5c_outputs(
        config,
        candidate_loci,
        phase5_locus,
        phase5b_inventory,
        integrity_output,
        locus_output,
        gene_output,
        report_output,
    )

    integrity = read_tsv(integrity_output, PHASE5C_INTEGRITY_COLUMNS)
    locus = read_tsv(locus_output, PHASE5C_LOCUS_QC_COLUMNS)[0]
    gene = read_tsv(gene_output, PHASE5C_GENE_QC_COLUMNS)[0]

    assert integrity[0]["candidate_context_usability"] == "CANDIDATE_LOCUS_CONTEXT_ONLY_WITH_PROCESSREPEATS_WARNING"
    assert any(row["file_role"] == "repeatmasker_tbl" and row["status"] == "MISSING" for row in integrity)
    assert locus["gff_out_overlap_concordance"] == "BOTH_GFF_AND_OUT_DIRECT_OVERLAP"
    assert locus["artifact_risk_modifier"] == "high_repeat_context_artifact_risk"
    assert gene["phase5c_qc_status"] == "REPEAT_CONTEXT_QC_SUPPORTED_BY_GFF_AND_OUT_ARTIFACT_ONLY"
    assert "artifact/context QC only" in gene["conservative_interpretation"]
