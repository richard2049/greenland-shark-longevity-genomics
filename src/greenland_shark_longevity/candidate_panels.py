"""Candidate panel validation."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .schemas import ALLOWED_MECHANISMS
from .utils import configure_logging, read_yaml, write_tsv

VALIDATION_COLUMNS = ["panel_id", "mechanism", "gene_symbol", "status", "message"]


def iter_candidates(panel_config: dict) -> list[dict[str, object]]:
    panels = panel_config.get("panels")
    if not isinstance(panels, dict):
        raise ValueError("Candidate panel YAML must contain a mapping named panels")
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    for panel_id, panel in panels.items():
        if not isinstance(panel, dict):
            raise ValueError(f"Panel {panel_id} must be a mapping")
        mechanism = panel.get("mechanism")
        candidates = panel.get("candidates")
        if mechanism not in ALLOWED_MECHANISMS:
            raise ValueError(f"Panel {panel_id} has unsupported mechanism {mechanism!r}")
        if not isinstance(candidates, list) or not candidates:
            raise ValueError(f"Panel {panel_id} must contain at least one candidate")
        for candidate in candidates:
            if not isinstance(candidate, dict):
                raise ValueError(f"Panel {panel_id} contains a non-mapping candidate")
            gene_symbol = str(candidate.get("gene_symbol", "")).strip()
            caveats = str(candidate.get("caveats", "")).strip()
            reference_ids = candidate.get("reference_ids")
            if not gene_symbol:
                raise ValueError(f"Panel {panel_id} contains a candidate without gene_symbol")
            key = (str(mechanism), gene_symbol.upper())
            if key in seen:
                raise ValueError(f"Duplicate candidate {gene_symbol} in mechanism {mechanism}")
            seen.add(key)
            if not caveats:
                raise ValueError(f"Candidate {gene_symbol} in {panel_id} must include caveats")
            if not isinstance(reference_ids, dict):
                raise ValueError(f"Candidate {gene_symbol} in {panel_id} must include reference_ids mapping")
            rows.append(
                {
                    "panel_id": panel_id,
                    "mechanism": mechanism,
                    "gene_symbol": gene_symbol,
                    "status": "OK",
                    "message": "Candidate schema valid.",
                }
            )
    return rows


def validate_panel_file(path: Path) -> list[dict[str, object]]:
    return iter_candidates(read_yaml(path))


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate candidate gene panels.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    configure_logging()
    logging.info("Validating candidate panels from %s", args.input)
    rows = validate_panel_file(args.input)
    write_tsv(args.output, rows, VALIDATION_COLUMNS)


if __name__ == "__main__":
    main()

