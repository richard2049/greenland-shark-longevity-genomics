from greenland_shark_longevity.hmmer_pfam import build_domain_integrity_rows, parse_domtblout
from greenland_shark_longevity.schemas import DOMAIN_CHECK_TARGET_COLUMNS
from greenland_shark_longevity.utils import write_tsv


def write_config(path):
    path.write_text(
        """
domain_validation:
  hmmer_executable: hmmscan
  pfam_hmm: data/raw/references/PFAM/Pfam-A.hmm
  threads: 1
  use_cut_ga: true
  max_i_evalue: 1.0e-5
  min_partial_hmm_coverage: 0.20
  min_full_hmm_coverage: 0.70
  recommended_command: hmmscan --cut_ga --domtblout out.domtblout Pfam-A.hmm candidates.faa
""",
        encoding="utf-8",
    )


def test_parse_domtblout_reads_hmmer_fields(tmp_path):
    domtblout = tmp_path / "pfam.domtblout"
    domtblout.write_text(
        "# HMMER domtblout fixture\n"
        "PF_repair PF00001.1 100 ERCC1|gene-1|prot1 - 315 1e-40 150.0 0.0 1 1 1e-42 1e-40 140.0 0.0 1 90 10 99 5 105 0.95 DNA repair domain\n",
        encoding="utf-8",
    )
    hits = parse_domtblout(domtblout)
    assert len(hits) == 1
    assert hits[0].accession_root == "PF00001"
    assert hits[0].query_name == "ERCC1|gene-1|prot1"
    assert round(hits[0].hmm_coverage, 2) == 0.90


def test_domain_integrity_classification_is_conservative(tmp_path):
    config = tmp_path / "config.yaml"
    write_config(config)
    targets = tmp_path / "domain_targets.tsv"
    write_tsv(
        targets,
        [
            {
                "mechanism": "DNA repair/genome stability",
                "gene_symbol": "ERCC1",
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "gene_id": "gene-1",
                "orthogroup_id": "OG1",
                "representative_protein_id": "prot1",
                "representative_fasta_id": "ERCC1|gene-1|prot1",
                "representative_length_aa": "315",
                "candidate_fasta": "candidates.faa",
                "domain_check_status": "READY_FOR_DOMAIN_SCAN",
                "recommended_method": "HMMER/Pfam",
                "required_validation": "REQUIRES_VALIDATION",
                "notes": "ready",
            },
            {
                "mechanism": "p53 pathway",
                "gene_symbol": "MDM2",
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "gene_id": "gene-2",
                "orthogroup_id": "OG2",
                "representative_protein_id": "prot2",
                "representative_fasta_id": "MDM2|gene-2|prot2",
                "representative_length_aa": "500",
                "candidate_fasta": "candidates.faa",
                "domain_check_status": "READY_FOR_DOMAIN_SCAN",
                "recommended_method": "HMMER/Pfam",
                "required_validation": "REQUIRES_VALIDATION",
                "notes": "ready",
            },
            {
                "mechanism": "Antioxidant response",
                "gene_symbol": "SOD2",
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "gene_id": "gene-3",
                "orthogroup_id": "OG3",
                "representative_protein_id": "prot3",
                "representative_fasta_id": "SOD2|gene-3|prot3",
                "representative_length_aa": "306",
                "candidate_fasta": "candidates.faa",
                "domain_check_status": "READY_FOR_DOMAIN_SCAN",
                "recommended_method": "HMMER/Pfam",
                "required_validation": "REQUIRES_VALIDATION",
                "notes": "ready",
            },
            {
                "mechanism": "Telomere-related biology",
                "gene_symbol": "TERT",
                "resource_id": "SMIC_TOKYO_GENOME_2025",
                "gene_id": "gene-4",
                "orthogroup_id": "OG4",
                "representative_protein_id": "NOT_ASSESSED",
                "representative_fasta_id": "NOT_ASSESSED",
                "representative_length_aa": "NOT_ASSESSED",
                "candidate_fasta": "candidates.faa",
                "domain_check_status": "SEQUENCE_UNAVAILABLE",
                "recommended_method": "HMMER/Pfam",
                "required_validation": "REQUIRES_VALIDATION",
                "notes": "missing",
            },
        ],
        DOMAIN_CHECK_TARGET_COLUMNS,
    )
    domtblout = tmp_path / "pfam.domtblout"
    domtblout.write_text(
        "PF_repair PF00001.1 100 ERCC1|gene-1|prot1 - 315 1e-40 150.0 0.0 1 1 1e-42 1e-40 140.0 0.0 1 90 10 99 5 105 0.95 DNA repair domain\n"
        "PF_ring PF00002.1 100 MDM2|gene-2|prot2 - 500 1e-20 90.0 0.0 1 1 1e-21 1e-20 80.0 0.0 1 30 250 279 245 284 0.80 Partial RING domain\n",
        encoding="utf-8",
    )

    rows = build_domain_integrity_rows(targets, domtblout, config)
    by_gene = {row["gene_symbol"]: row for row in rows}
    assert by_gene["ERCC1"]["domain_validation_status"] == "DOMAIN_SUPPORTED"
    assert by_gene["MDM2"]["domain_validation_status"] == "PARTIAL_DOMAIN"
    assert by_gene["SOD2"]["domain_validation_status"] == "NO_EXPECTED_DOMAIN_DETECTED"
    assert by_gene["TERT"]["domain_validation_status"] == "NOT_ASSESSED"
    assert "not gene absence" in by_gene["SOD2"]["classification_rule"]
