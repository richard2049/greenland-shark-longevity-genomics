from pathlib import Path

import pytest

from greenland_shark_longevity.metadata import build_manifests, validate_resources
from greenland_shark_longevity.schemas import DATA_MANIFEST_COLUMNS
from greenland_shark_longevity.utils import read_tsv, read_yaml


def test_metadata_manifest_generation(tmp_path):
    config = read_yaml(Path("config/config.yaml"))
    paths = build_manifests(config, tmp_path)
    data_rows = read_tsv(paths["data_manifest"], DATA_MANIFEST_COLUMNS)
    status_rows = read_tsv(paths["resource_status"])
    assert {row["resource_id"] for row in data_rows} >= {
        "SMIC_FLI_GENOME_2025",
        "SMIC_TOKYO_GENOME_2025",
        "SMIC_RETINA_PRJNA1246101_2026",
    }
    assert {row["status"] for row in status_rows}.issubset({"REGISTERED", "LOCAL_AVAILABLE", "MISSING_LOCAL"})


def test_invalid_resource_status_fails():
    resource = {
        "resource_id": "bad",
        "species": "Somniosus microcephalus",
        "accession": "TODO",
        "source": "TODO",
        "resource_type": "genome",
        "status": "UNKNOWN",
        "size": "TODO",
        "local_path": "TODO",
        "source_url": "TODO",
        "publication_or_resource": "TODO",
        "usage_notes": "TODO",
        "expected_files": "genome_fasta",
    }
    with pytest.raises(ValueError, match="invalid status"):
        validate_resources([resource])
