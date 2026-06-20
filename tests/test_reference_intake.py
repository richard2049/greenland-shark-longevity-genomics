from pathlib import Path

from greenland_shark_longevity.reference_intake import (
    configured_reference_files,
    parse_assembly_stats,
    parse_feature_count,
    parse_gff_gene_counts,
)
from greenland_shark_longevity.utils import read_yaml


def test_parse_assembly_stats_prefers_primary_scope(tmp_path):
    stats = tmp_path / "assembly_stats.txt"
    stats.write_text(
        "\n".join(
            [
                "# Assembly name:  ExampleAssembly",
                "# BioProject:     PRJNA000000",
                "# BioSample:      SAMN000000",
                "# Assembly level: Chromosome",
                "# Date:           2026-05-27",
                "# unit-name\tmolecule-name\tmolecule-type/loc\tsequence-type\tstatistic\tvalue",
                "all\tall\tall\tall\ttotal-length\t200",
                "all\tall\tall\tall\tscaffold-count\t4",
                "Primary Assembly\tall\tall\tall\ttotal-length\t180",
                "Primary Assembly\tall\tall\tall\tungapped-length\t170",
                "Primary Assembly\tall\tall\tall\tscaffold-count\t3",
                "Primary Assembly\tall\tall\tall\tscaffold-N50\t90",
                "Primary Assembly\tall\tall\tall\tgc-perc\t48",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    metadata, parsed, scope = parse_assembly_stats(stats)
    assert scope == "Primary Assembly"
    assert metadata["assembly name"] == "ExampleAssembly"
    assert parsed["total-length"] == "180"
    assert parsed["scaffold-N50"] == "90"


def test_parse_feature_count_counts_genes(tmp_path):
    feature_count = tmp_path / "feature_count.txt"
    feature_count.write_text(
        "\n".join(
            [
                "# Feature\tClass\tFull Assembly\tAssembly-unit accession\tAssembly-unit name\tUnique Ids\tPlacements",
                "gene\tprotein_coding\tGCA_000000000.1\tGCA_000000001.1\tPrimary Assembly\t7\t7",
                "gene\tprotein_coding\tGCA_000000000.1\t\tall\t7\t7",
                "gene\tncRNA\tGCA_000000000.1\t\tall\t3\t3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protein_coding, total_gene = parse_feature_count(feature_count)
    assert protein_coding == 7
    assert total_gene == 10


def test_parse_gff_gene_counts(tmp_path):
    gff = tmp_path / "annotation.gff"
    gff.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "scaffold_1\tGnomon\tgene\t10\t100\t.\t+\t.\tID=gene-a;gene_biotype=protein_coding",
                "scaffold_1\tGnomon\tmRNA\t10\t100\t.\t+\t.\tID=rna-a;Parent=gene-a",
                "scaffold_1\tGnomon\tgene\t200\t300\t.\t-\t.\tID=gene-b;gene_biotype=lncRNA",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protein_coding, total_gene = parse_gff_gene_counts(gff)
    assert protein_coding == 1
    assert total_gene == 2


def test_comparator_refseq_protein_sources_are_registered():
    config = read_yaml(Path("config/config.yaml"))
    rows = configured_reference_files(config)
    by_resource = {row["resource_id"]: row for row in rows if row["file_role"] == "protein_fasta"}
    for resource_id in ["CMIL_REFSEQ_2021", "SCAN_REFSEQ_2026", "ARAD_REFSEQ_2026", "RTYP_REFSEQ_2022"]:
        assert by_resource[resource_id]["url"].endswith("_protein.faa.gz")
        assert by_resource[resource_id]["selected_for_download"] == "True"
