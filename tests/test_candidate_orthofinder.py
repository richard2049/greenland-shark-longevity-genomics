from pathlib import Path

from greenland_shark_longevity.candidate_orthofinder import (
    integrate_real_candidates,
    parse_gff_protein_annotations,
)
from greenland_shark_longevity.schemas import (
    ORTHOFINDER_INPUT_MANIFEST_COLUMNS,
    ORTHOGROUP_GENE_COUNT_LONG_COLUMNS,
    REFERENCE_GENE_COORDINATE_COLUMNS,
)
from greenland_shark_longevity.utils import write_tsv


def test_parse_gff_protein_annotations_collapses_cds_exons(tmp_path):
    gff = tmp_path / "annotation.gff"
    gff.write_text(
        "scaffold_1\tGnomon\tCDS\t1\t30\t.\t+\t0\t"
        "ID=cds-gs_001-P1;protein_id=gnl|WGS:ZZZZ|gs_001-P1;gene=FTH1;locus_tag=gs_001;product=ferritin heavy chain\n"
        "scaffold_1\tGnomon\tCDS\t40\t60\t.\t+\t0\t"
        "ID=cds-gs_001-P1;protein_id=gnl|WGS:ZZZZ|gs_001-P1;gene=FTH1;locus_tag=gs_001;product=ferritin heavy chain\n",
        encoding="utf-8",
    )
    annotations = parse_gff_protein_annotations(gff)
    assert len(annotations) == 1
    assert annotations[0].gene_id == "gene-gs_001"
    assert annotations[0].gene_symbol == "FTH1"


def test_real_candidate_mapping_uses_exact_symbols_and_collapses_isoforms(tmp_path):
    panels = tmp_path / "candidate_panels.yaml"
    panels.write_text(
        """
metadata:
  schema_version: "test"
panels:
  ferroptosis:
    mechanism: Ferroptosis, iron handling, and oxidative stress
    candidates:
      - gene_symbol: FTH1
        synonyms: [ferritin heavy chain]
        reference_ids: {human_gene: TODO}
        caveats: Test caveat.
      - gene_symbol: FTH1B
        synonyms: [FTH1b, ferritin heavy chain 1b]
        reference_ids: {human_gene: NOT_ASSESSED}
        caveats: Test caveat.
""",
        encoding="utf-8",
    )
    gff = tmp_path / "annotation.gff"
    gff.write_text(
        "\n".join(
            [
                "scaffold_1\tGnomon\tCDS\t1\t30\t.\t+\t0\tID=cds-gs_001-P1;protein_id=gnl|WGS:ZZZZ|gs_001-P1;gene=FTH1;locus_tag=gs_001;product=ferritin heavy chain",
                "scaffold_1\tGnomon\tCDS\t40\t70\t.\t+\t0\tID=cds-gs_001-P2;protein_id=gnl|WGS:ZZZZ|gs_001-P2;gene=FTH1;locus_tag=gs_001;product=ferritin heavy chain isoform 2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    result_root = tmp_path / "Results_Test"
    orthogroups = result_root / "Orthogroups"
    orthogroups.mkdir(parents=True)
    (orthogroups / "Orthogroups.tsv").write_text(
        "Orthogroup\tcmil__CMIL_REFSEQ_2021\tsmic__SMIC_TOKYO_GENOME_2025\n"
        "OG0000001\tXP_1.1\tgnl|WGS_ZZZZ|gs_001-P1, gnl|WGS_ZZZZ|gs_001-P2\n",
        encoding="utf-8",
    )
    gene_counts = tmp_path / "orthogroup_gene_counts_long.tsv"
    write_tsv(
        gene_counts,
        [
            {
                "orthogroup_id": "OG0000001",
                "species_id": "cmil",
                "resource_id": "CMIL_REFSEQ_2021",
                "orthofinder_species_column": "cmil__CMIL_REFSEQ_2021",
                "copy_count": "1",
                "orthofinder_results_dir": str(result_root),
            },
            {
                "orthogroup_id": "OG0000001",
                "species_id": "smic",
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "orthofinder_species_column": "smic__SMIC_TOKYO_GENOME_2025",
                "copy_count": "2",
                "orthofinder_results_dir": str(result_root),
            },
        ],
        ORTHOGROUP_GENE_COUNT_LONG_COLUMNS,
    )
    coordinates = tmp_path / "coordinates.tsv"
    write_tsv(
        coordinates,
        [
            {
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "accession": "GCA_TEST",
                "gene_id": "gene-gs_001",
                "gene_symbol": "FTH1",
                "seqid": "scaffold_1",
                "start": "1",
                "end": "100",
                "strand": "+",
                "source_file": str(gff),
                "parse_status": "PARSED_FROM_GFF",
                "notes": "test",
            }
        ],
        REFERENCE_GENE_COORDINATE_COLUMNS,
    )
    manifest = tmp_path / "manifest.tsv"
    write_tsv(
        manifest,
        [
            {
                "species_id": "smic",
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "source_protein_fasta": "smic.faa",
                "orthofinder_input_path": "smic__SMIC_TOKYO_GENOME_2025.faa",
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
                "resource_id": "CMIL_REFSEQ_2021",
                "source_protein_fasta": "cmil.faa",
                "orthofinder_input_path": "cmil__CMIL_REFSEQ_2021.faa",
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

    copy_rows, audit_rows = integrate_real_candidates(panels, gff, gene_counts, coordinates, manifest)
    copy_by_gene = {row["gene_symbol"]: row for row in copy_rows}
    audit_by_gene = {row["gene_symbol"]: row for row in audit_rows}

    assert copy_by_gene["FTH1"]["copy_count"] == "1"
    assert copy_by_gene["FTH1"]["orthogroup_target_protein_count"] == "2"
    assert audit_by_gene["FTH1"]["isoform_risk"] == "yes"
    assert audit_by_gene["FTH1"]["orthology_support"] == "yes"
    assert copy_by_gene["FTH1B"]["mapping_status"] == "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH"
    assert copy_by_gene["FTH1B"]["copy_count"] == "0"
