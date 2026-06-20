from pathlib import Path

from greenland_shark_longevity.phase5_repeat_context import (
    build_candidate_windows,
    build_gene_context_rows,
    build_locus_context_rows,
    discover_repeat_resources,
    load_repeat_features,
    load_assembly_seqid_aliases,
    parse_gff_repeat_features,
)
from greenland_shark_longevity.schemas import LOCAL_FILE_INVENTORY_COLUMNS, PHASE4E_LOCUS_HARDENING_COLUMNS
from greenland_shark_longevity.utils import write_tsv


def write_locus_table(path: Path) -> None:
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


def test_discovery_records_no_repeat_annotation_when_figshare_inventory_has_only_gene_files(tmp_path: Path):
    gff = tmp_path / "complete.genomic.gff"
    gff.write_text(
        "##gff-version 3\n"
        "scaffold_33\tGnomon\tgene\t10\t100\t.\t+\t.\tID=gene-a;gene_biotype=protein_coding\n",
        encoding="utf-8",
    )
    inventory = tmp_path / "inventory.tsv"
    write_tsv(
        inventory,
        [
            {
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "source_id": "PNAS_YANG_2026_FIGSHARE_ANNOTATION",
                "file_name": "complete.genomic.gff",
                "file_role_candidate": "annotation_gff",
                "detected_format": "gff3",
                "local_path": str(gff),
                "byte_size": str(gff.stat().st_size),
                "sequence_count": "NOT_APPLICABLE",
                "notes": "gene annotation",
            }
        ],
        LOCAL_FILE_INVENTORY_COLUMNS,
    )

    rows = discover_repeat_resources(inventory)

    assert any(row["repeat_annotation_status"] == "NO_REPEAT_FEATURES_DETECTED" for row in rows)
    assert rows[-1]["repeat_annotation_status"] == "NO_REPEAT_ANNOTATION_FOUND"


def test_gff_repeat_parser_recognizes_repeat_region_and_mobile_element(tmp_path: Path):
    gff = tmp_path / "repeats.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "scaffold_33\tRepeatMasker\trepeat_region\t100\t150\t.\t+\t.\tID=r1;repeat_class=LINE;repeat_family=L2;Name=L2-1\n"
        "scaffold_33\tRepeatMasker\tmobile_element\t180\t220\t.\t-\t.\tID=r2;mobile_element_type=LTR_retrotransposon;Name=ERV1\n",
        encoding="utf-8",
    )

    features = parse_gff_repeat_features(gff)

    assert len(features) == 2
    assert features[0].repeat_class == "LINE"
    assert features[0].repeat_family == "L2"
    assert features[1].repeat_class == "LTR_retrotransposon"


def test_gff_repeat_parser_recognizes_repeatmasker_gff_similarity(tmp_path: Path):
    gff = tmp_path / "smic_tokyo.repeatmasker.out.gff"
    gff.write_text(
        "##gff-version 3\n"
        "scaffold_33\tRepeatMasker\tsimilarity\t100\t150\t22.0\t+\t.\tTarget=Motif:LINE_L2 1 51\n",
        encoding="utf-8",
    )

    features = parse_gff_repeat_features(gff, candidate_by_name=True)

    assert len(features) == 1
    assert features[0].seqid == "scaffold_33"
    assert features[0].repeat_name == "LINE_L2"


def test_locus_intersection_uses_not_assessed_when_no_repeat_annotation_is_available(tmp_path: Path):
    loci = tmp_path / "loci.tsv"
    write_locus_table(loci)
    resource_rows = [
        {
            "repeat_annotation_status": "NO_REPEAT_ANNOTATION_FOUND",
        }
    ]

    context = build_locus_context_rows(loci, [], resource_rows, window_bp=10)

    assert context[0]["repeat_annotation_status"] == "NO_REPEAT_ANNOTATION_AVAILABLE"
    assert context[0]["repeat_overlap_count"] == "NOT_ASSESSED"
    assert context[0]["artifact_context_status"] == "REPEAT_CONTEXT_NOT_ASSESSABLE_CURRENT_LOCAL_FILES"


def test_locus_intersection_records_overlap_when_repeat_annotations_exist(tmp_path: Path):
    loci = tmp_path / "loci.tsv"
    write_locus_table(loci)
    repeat_gff = tmp_path / "repeat_annotation.gff3"
    repeat_gff.write_text(
        "##gff-version 3\n"
        "scaffold_33\tRepeatMasker\trepeat_region\t150\t250\t.\t+\t.\tID=r1;repeat_class=LINE;repeat_family=L2;Name=L2-1\n",
        encoding="utf-8",
    )
    inventory = tmp_path / "inventory.tsv"
    write_tsv(
        inventory,
        [
            {
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "source_id": "TEST",
                "file_name": "repeat_annotation.gff3",
                "file_role_candidate": "repeat_annotation",
                "detected_format": "gff3",
                "local_path": str(repeat_gff),
                "byte_size": str(repeat_gff.stat().st_size),
                "sequence_count": "NOT_APPLICABLE",
                "notes": "test repeat annotation",
            }
        ],
        LOCAL_FILE_INVENTORY_COLUMNS,
    )

    resource_rows = discover_repeat_resources(inventory)
    features = load_repeat_features(resource_rows)
    context = build_locus_context_rows(loci, features, resource_rows, window_bp=10)
    summary = build_gene_context_rows(context)

    assert context[0]["repeat_annotation_status"] == "REPEAT_ANNOTATION_AVAILABLE"
    assert context[0]["repeat_overlap_count"] == "1"
    assert context[0]["repeat_overlap_bp"] == "51"
    assert context[0]["repeat_classes"] == "LINE"
    assert summary[0]["loci_with_repeat_overlap"] == "1"


def test_repeatmasker_accessions_are_normalized_and_filtered_to_candidate_windows(tmp_path: Path):
    loci = tmp_path / "loci.tsv"
    write_locus_table(loci)
    assembly_report = tmp_path / "assembly_report.txt"
    assembly_report.write_text(
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
        "scaffold_33\tunplaced-scaffold\tna\tna\tJBLTJD010000033.1\t<>\tna\tPrimary Assembly\t1000000\tna\n",
        encoding="utf-8",
    )
    repeat_gff = tmp_path / "smic_tokyo.repeatmasker.out.gff"
    repeat_gff.write_text(
        "##gff-version 3\n"
        "JBLTJD010000033.1\tRepeatMasker\tdispersed_repeat\t120\t180\t1.0\t+\t.\tID=r1;Target \"Motif:rnd-1_family-1\" 1 61\n"
        "JBLTJD010000033.1\tRepeatMasker\tdispersed_repeat\t500000\t500100\t1.0\t+\t.\tID=r2;Target \"Motif:rnd-1_family-2\" 1 101\n",
        encoding="utf-8",
    )
    inventory = tmp_path / "inventory.tsv"
    write_tsv(
        inventory,
        [
            {
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "source_id": "TEST",
                "file_name": "smic_tokyo.repeatmasker.out.gff",
                "file_role_candidate": "repeat_annotation",
                "detected_format": "gff3",
                "local_path": str(repeat_gff),
                "byte_size": str(repeat_gff.stat().st_size),
                "sequence_count": "NOT_APPLICABLE",
                "notes": "test repeat annotation",
            }
        ],
        LOCAL_FILE_INVENTORY_COLUMNS,
    )

    aliases = load_assembly_seqid_aliases(assembly_report)
    windows = build_candidate_windows(loci, window_bp=10, seqid_aliases=aliases)
    resource_rows = discover_repeat_resources(inventory, seqid_aliases=aliases, target_windows=windows)
    features = load_repeat_features(resource_rows, seqid_aliases=aliases, target_windows=windows)
    context = build_locus_context_rows(loci, features, resource_rows, window_bp=10)

    assert len(features) == 1
    assert features[0].seqid == "scaffold_33"
    assert "normalized to scaffold_33" in features[0].notes
    assert context[0]["repeat_overlap_count"] == "1"
    assert context[0]["repeat_overlap_bp"] == "61"
