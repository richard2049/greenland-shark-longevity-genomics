from pathlib import Path

from greenland_shark_longevity.phase3_inputs import (
    build_orthofinder_input_manifest,
    configured_source_checks,
    parse_gff_gene_coordinates,
    stage_orthofinder_input,
)
from greenland_shark_longevity.schemas import REFERENCE_FILE_COLUMNS, REFERENCE_PROTEIN_QC_COLUMNS
from greenland_shark_longevity.utils import read_tsv, read_yaml


def test_configured_source_checks_include_pnas_and_ncbi_searches():
    config = read_yaml(Path("config/config.yaml"))
    rows = configured_source_checks(config)
    by_id = {row["source_id"]: row for row in rows}
    assert "PNAS_YANG_2026_FINAL_ARTICLE" in by_id
    assert "PNAS_YANG_2026_FIGSHARE_ANNOTATION" in by_id
    assert "NCBI_DATASETS_GCA_056099535_1" in by_id
    assert by_id["PNAS_YANG_2026_FINAL_ARTICLE"]["selected_for_download"] == "False"
    assert "protein FASTA" in by_id["PNAS_YANG_2026_FINAL_ARTICLE"]["searched_for"]
    assert by_id["PNAS_YANG_2026_FIGSHARE_ANNOTATION"]["status"] == "LOCAL_FILES_MANUALLY_DOWNLOADED_INVENTORIED"


def test_orthofinder_manifest_excludes_non_genome_wide_cds_subset():
    config = read_yaml(Path("config/config.yaml"))
    inventory = read_tsv(Path("data/metadata/reference_file_inventory.tsv"), REFERENCE_FILE_COLUMNS)
    protein_qc = read_tsv(Path("results/qc/reference_protein_qc.tsv"), REFERENCE_PROTEIN_QC_COLUMNS)
    rows = build_orthofinder_input_manifest(config, inventory, protein_qc)
    by_resource = {row["resource_id"]: row for row in rows}
    assert by_resource["SMIC_RETINA_PRJNA1246101_2026"]["orthofinder_ready"] == "False"
    assert by_resource["SMIC_RETINA_PRJNA1246101_2026"]["blocker"] == "PROTEIN_FASTA_IS_NOT_GENOME_WIDE"
    assert by_resource["SMIC_FLI_GENOME_2025"]["blocker"] == "NO_GENOME_WIDE_PROTEIN_FASTA_REGISTERED"
    assert by_resource["SMIC_TOKYO_GENOME_2025"]["orthofinder_ready"] == "True"
    assert by_resource["SMIC_TOKYO_GENOME_2025"]["protein_set_scope"] == "GENOME_WIDE_PROTEIN_FASTA_CANDIDATE"


def test_parse_gff_gene_coordinates(tmp_path):
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "scaf1\tmaker\tgene\t10\t50\t.\t+\t.\tID=gene1;Name=TP53",
                "scaf1\tmaker\tmRNA\t10\t50\t.\t+\t.\tID=tx1;Parent=gene1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    rows = parse_gff_gene_coordinates(gff, "RESOURCE", "GCA_000000000.1")
    assert len(rows) == 1
    assert rows[0]["gene_id"] == "gene1"
    assert rows[0]["gene_symbol"] == "TP53"
    assert rows[0]["parse_status"] == "PARSED_FROM_GFF"


def test_stage_orthofinder_input_decompresses_gzip(tmp_path):
    import gzip

    source = tmp_path / "proteins.faa.gz"
    target = tmp_path / "orthofinder" / "species.faa"
    with gzip.open(source, "wb") as handle:
        handle.write(b">p1\nMA\n")
    stage_orthofinder_input(source, target)
    assert target.read_text(encoding="utf-8") == ">p1\nMA\n"
