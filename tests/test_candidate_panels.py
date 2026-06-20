from pathlib import Path

import pytest
import yaml

from greenland_shark_longevity.candidate_panels import iter_candidates, validate_panel_file
from greenland_shark_longevity.schemas import ALLOWED_MECHANISMS


def test_candidate_panel_fixture_validates():
    rows = validate_panel_file(Path("config/candidate_panels.yaml"))
    assert rows
    assert {row["status"] for row in rows} == {"OK"}
    assert {row["mechanism"] for row in rows}.issubset(ALLOWED_MECHANISMS)


def test_candidate_panel_requires_caveats(tmp_path):
    panel = {
        "panels": {
            "bad": {
                "mechanism": "DNA repair/genome stability",
                "candidates": [
                    {"gene_symbol": "TP53", "synonyms": [], "reference_ids": {"human_gene": "TODO"}}
                ],
            }
        }
    }
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(panel), encoding="utf-8")
    with pytest.raises(ValueError, match="caveats"):
        validate_panel_file(path)


def test_duplicate_gene_within_mechanism_fails():
    panel = {
        "panels": {
            "bad": {
                "mechanism": "p53 pathway",
                "candidates": [
                    {
                        "gene_symbol": "TP53",
                        "synonyms": [],
                        "reference_ids": {"human_gene": "TODO"},
                        "caveats": "first",
                    },
                    {
                        "gene_symbol": "tp53",
                        "synonyms": [],
                        "reference_ids": {"human_gene": "TODO"},
                        "caveats": "second",
                    },
                ],
            }
        }
    }
    with pytest.raises(ValueError, match="Duplicate candidate"):
        iter_candidates(panel)

