from greenland_shark_longevity.candidate_validation import (
    build_isoform_and_domain_rows,
    build_rescue_rows,
    parse_fasta,
)
from greenland_shark_longevity.schemas import COPY_NUMBER_COLUMNS, DUPLICATION_AUDIT_COLUMNS
from greenland_shark_longevity.utils import write_tsv


def test_candidate_validation_selects_longest_representative_isoform(tmp_path):
    panels = tmp_path / "candidate_panels.yaml"
    panels.write_text(
        """
metadata:
  schema_version: "test"
panels:
  dna_repair:
    mechanism: DNA repair/genome stability
    candidates:
      - gene_symbol: ERCC1
        synonyms: []
        reference_ids: {human_gene: TODO}
        caveats: Test caveat.
""",
        encoding="utf-8",
    )
    copy_number = tmp_path / "copy.tsv"
    write_tsv(
        copy_number,
        [
            {
                "mechanism": "DNA repair/genome stability",
                "gene_symbol": "ERCC1",
                "species_id": "smic",
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "orthogroup_id": "OG1",
                "copy_count": "1",
                "gene_ids": "gene-gs_001",
                "protein_ids": "gnl|WGS_ZZZZ|gs_001-P1,gnl|WGS_ZZZZ|gs_001-P2",
                "orthogroup_target_protein_count": "2",
                "orthogroup_species_counts": "OG1:smic=2",
                "annotation_match_level": "exact",
                "mapping_status": "ANNOTATION_SYMBOL_MATCH_ORTHOGROUP_MAPPED",
                "orthofinder_results_dir": "results",
                "resources_supporting": "SMIC_TOKYO_GENOME_2025",
                "demo_only": "False",
            }
        ],
        COPY_NUMBER_COLUMNS,
    )
    audit = tmp_path / "audit.tsv"
    write_tsv(
        audit,
        [
            {
                "mechanism": "DNA repair/genome stability",
                "gene_symbol": "ERCC1",
                "orthogroup_id": "OG1",
                "copy_count": "1",
                "copy_ids": "gene-gs_001",
                "protein_ids": "gnl|WGS_ZZZZ|gs_001-P1,gnl|WGS_ZZZZ|gs_001-P2",
                "orthology_support": "yes",
                "domain_integrity": "NOT_ASSESSED",
                "separable_loci": "no",
                "coordinate_support": "yes",
                "cross_resource_support": "no",
                "isoform_risk": "yes",
                "fragmentation_risk": "no",
                "expression_support": "NOT_ASSESSED",
                "annotation_match_level": "exact",
                "orthogroup_species_counts": "OG1:smic=2",
                "coordinate_summary": "gene-gs_001:scaffold_1:1-100:+",
                "mapping_status": "ANNOTATION_SYMBOL_MATCH_ORTHOGROUP_MAPPED",
                "orthofinder_results_dir": "results",
                "resources_supporting": "SMIC_TOKYO_GENOME_2025",
                "artifact_risk": "moderate",
                "required_validation": "REQUIRES_VALIDATION",
                "demo_only": "False",
            }
        ],
        DUPLICATION_AUDIT_COLUMNS,
    )
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(
        ">gnl|WGS:ZZZZ|gs_001-P1\nMAAA\n"
        ">gnl|WGS:ZZZZ|gs_001-P2\nMAAAAAA\n",
        encoding="utf-8",
    )
    representative_fasta = tmp_path / "representatives.faa"
    isoform_rows, domain_rows, sequences = build_isoform_and_domain_rows(
        panels, copy_number, audit, proteins, representative_fasta
    )
    assert isoform_rows[0]["isoform_count"] == "2"
    assert isoform_rows[0]["representative_protein_id"] == "gnl|WGS_ZZZZ|gs_001-P2"
    assert isoform_rows[0]["representative_length_aa"] == "7"
    assert domain_rows[0]["domain_check_status"] == "READY_FOR_DOMAIN_SCAN"
    assert list(sequences.values()) == ["MAAAAAA"]


def test_rescue_rows_keep_unresolved_genes_as_annotation_uncertainty(tmp_path):
    panels = tmp_path / "candidate_panels.yaml"
    panels.write_text(
        """
metadata:
  schema_version: "test"
panels:
  p53:
    mechanism: p53 pathway
    candidates:
      - gene_symbol: TP53
        synonyms: [p53]
        reference_ids: {human_gene: TODO}
        caveats: Test caveat.
""",
        encoding="utf-8",
    )
    copy_number = tmp_path / "copy.tsv"
    write_tsv(
        copy_number,
        [
            {
                "mechanism": "p53 pathway",
                "gene_symbol": "TP53",
                "species_id": "smic",
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "orthogroup_id": "NOT_ASSESSED",
                "copy_count": "0",
                "gene_ids": "NOT_ASSESSED",
                "protein_ids": "NOT_ASSESSED",
                "orthogroup_target_protein_count": "NOT_ASSESSED",
                "orthogroup_species_counts": "NOT_ASSESSED",
                "annotation_match_level": "NOT_ASSESSED",
                "mapping_status": "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH",
                "orthofinder_results_dir": "results",
                "resources_supporting": "NOT_ASSESSED",
                "demo_only": "False",
            }
        ],
        COPY_NUMBER_COLUMNS,
    )
    rows = build_rescue_rows(panels, copy_number, ["TP53"])
    assert rows[0]["gene_symbol"] == "TP53"
    assert rows[0]["mapping_status"] == "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH"
    assert "do not infer absence" in rows[0]["blockers"]


def test_parse_fasta_uses_first_token_as_identifier(tmp_path):
    fasta = tmp_path / "x.faa"
    fasta.write_text(">protein1 description\nMA\n", encoding="utf-8")
    assert parse_fasta(fasta) == {"protein1": "MA"}
