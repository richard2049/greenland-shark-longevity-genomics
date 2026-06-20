from pathlib import Path

from greenland_shark_longevity.annotation_rescue import (
    FastaRecord,
    record_matches,
    select_query_records,
    summary_rows,
)


def test_phase4_query_patterns_are_gene_specific():
    h1 = FastaRecord("q1", "histone H1.0-like [test]", "MAAA", "x.faa", "scan", "SCAN")
    h101 = FastaRecord("q2", "histone H1.01-like [test]", "MAAA", "x.faa", "scan", "SCAN")
    rad51 = FastaRecord("q3", "DNA repair protein RAD51 homolog 1 [test]", "MAAA", "x.faa", "scan", "SCAN")
    rad51_paralog = FastaRecord("q4", "DNA repair protein RAD51 homolog 3 [test]", "MAAA", "x.faa", "scan", "SCAN")
    assert record_matches(h1, "H1F0")
    assert not record_matches(h101, "H1F0")
    assert record_matches(rad51, "RAD51")
    assert not record_matches(rad51_paralog, "RAD51")


def test_select_query_records_limits_per_species():
    targets = [
        {
            "gene_symbol": "FTH1B",
            "mechanism": "Ferroptosis, iron handling, and oxidative stress",
            "mapping_status": "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH",
        }
    ]
    proteomes = {
        "smic": [FastaRecord("target", "ferritin heavy chain B-like", "MAAA", "smic.faa", "smic", "SMIC")],
        "scan": [
            FastaRecord("q1", "ferritin heavy chain B-like", "MAAA", "scan.faa", "scan", "SCAN"),
            FastaRecord("q2", "ferritin heavy chain B-like", "MAAAAA", "scan.faa", "scan", "SCAN"),
            FastaRecord("q3", "ferritin heavy chain B-like", "MA", "scan.faa", "scan", "SCAN"),
        ],
    }
    rows, selected = select_query_records(targets, proteomes, max_per_species=2)
    assert len(selected["FTH1B"]) == 2
    assert selected["FTH1B"][0].protein_id == "q2"
    assert [row["selection_status"] for row in rows].count("MATCHED_BUT_EXCLUDED_BY_QUERY_CAP") == 1


def test_phase4_summary_never_reports_absence_for_no_rescue():
    targets = [
        {
            "gene_symbol": "TP53",
            "mechanism": "p53 pathway",
            "mapping_status": "ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH",
        }
    ]
    rows = summary_rows(
        targets=targets,
        selected_queries={"TP53": []},
        forward_hits=[],
        reciprocal_supported={},
        coordinates={},
        genome_fasta=Path("missing.fa"),
    )
    assert rows[0]["protein_rescue_status"] == "QUERY_NOT_AVAILABLE"
    assert rows[0]["genome_alignment_status"] == "GENOME_FASTA_MISSING_NOT_RUN"
    assert "absent" in rows[0]["claim_language_guardrail"]
    assert rows[0]["classification"] == "NOT_ASSESSED"
