"""Metadata and provenance manifest generation."""

from __future__ import annotations

import argparse
import logging
from datetime import datetime, timezone
from pathlib import Path

from .schemas import (
    ALLOWED_STATUSES,
    DATA_MANIFEST_COLUMNS,
    RESOURCE_STATUS_COLUMNS,
    SPECIES_MANIFEST_COLUMNS,
)
from .utils import configure_logging, read_yaml, write_tsv


def _normalize_list(values: object) -> str:
    if values is None:
        return ""
    if isinstance(values, list):
        return ";".join(str(value) for value in values)
    return str(values)


def validate_resources(resources: list[dict]) -> None:
    required = set(DATA_MANIFEST_COLUMNS + ["expected_files"])
    for resource in resources:
        missing = sorted(required - set(resource))
        if missing:
            raise ValueError(f"Resource {resource.get('resource_id', '<unknown>')} missing fields: {', '.join(missing)}")
        if resource["status"] not in ALLOWED_STATUSES:
            raise ValueError(f"Resource {resource['resource_id']} has invalid status {resource['status']}")


def validate_species(species: list[dict]) -> None:
    required = set(SPECIES_MANIFEST_COLUMNS)
    for entry in species:
        missing = sorted(required - set(entry))
        if missing:
            raise ValueError(f"Species {entry.get('species_id', '<unknown>')} missing fields: {', '.join(missing)}")
        if entry["status"] not in ALLOWED_STATUSES:
            raise ValueError(f"Species {entry['species_id']} has invalid status {entry['status']}")


def build_manifests(config: dict, output_dir: Path) -> dict[str, Path]:
    resources = config.get("resources", [])
    species = config.get("species", [])
    if not isinstance(resources, list) or not isinstance(species, list):
        raise ValueError("Config must contain list-valued resources and species sections")

    validate_resources(resources)
    validate_species(species)

    output_dir.mkdir(parents=True, exist_ok=True)
    data_rows = [{column: _normalize_list(resource.get(column, "")) for column in DATA_MANIFEST_COLUMNS} for resource in resources]
    species_rows = [{column: _normalize_list(entry.get(column, "")) for column in SPECIES_MANIFEST_COLUMNS} for entry in species]
    checked_at = datetime.now(timezone.utc).isoformat()
    status_rows = []
    for resource in resources:
        local_path = str(resource.get("local_path", "TODO"))
        if local_path in {"", "TODO", "NOT_ASSESSED"}:
            status = "REGISTERED"
            message = "Registered public resource; no local file expected in MVP."
        else:
            exists = Path(local_path).exists()
            status = "LOCAL_AVAILABLE" if exists else "MISSING_LOCAL"
            message = "Expected local file found." if exists else "Expected local file is missing."
        status_rows.append(
            {
                "resource_id": resource["resource_id"],
                "expected_files": _normalize_list(resource["expected_files"]),
                "local_path": local_path,
                "status": status,
                "checked_at": checked_at,
                "message": message,
            }
        )

    paths = {
        "data_manifest": output_dir / "data_manifest.tsv",
        "species_manifest": output_dir / "species_manifest.tsv",
        "resource_status": output_dir / "resource_status.tsv",
    }
    write_tsv(paths["data_manifest"], data_rows, DATA_MANIFEST_COLUMNS)
    write_tsv(paths["species_manifest"], species_rows, SPECIES_MANIFEST_COLUMNS)
    write_tsv(paths["resource_status"], status_rows, RESOURCE_STATUS_COLUMNS)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate metadata/provenance manifests.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    configure_logging()
    logging.info("Generating metadata manifests from %s", args.config)
    build_manifests(read_yaml(args.config), args.output_dir)


if __name__ == "__main__":
    main()

