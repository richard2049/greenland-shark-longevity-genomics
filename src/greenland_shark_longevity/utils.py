"""Small IO and validation helpers."""

from __future__ import annotations

import csv
import gzip
import logging
from pathlib import Path
from typing import Iterable, TextIO

import yaml

NOT_ASSESSED = "NOT_ASSESSED"
MISSING_VALUE_MARKERS = {"", "TODO", NOT_ASSESSED}


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def open_text(path: Path, errors: str = "replace") -> TextIO:
    """Open plain or gzipped text files with UTF-8 decoding."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors=errors)
    return path.open("r", encoding="utf-8", errors=errors)


def clean_text(value: object | None) -> str:
    if value is None:
        return ""
    return str(value).strip()


def read_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected YAML mapping in {path}")
    return data


def write_tsv(path: Path, rows: Iterable[dict], columns: list[str]) -> None:
    ensure_parent(path)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path, required_columns: list[str] | None = None) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header")
        if required_columns:
            missing = [column for column in required_columns if column not in reader.fieldnames]
            if missing:
                raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
        return [dict(row) for row in reader]


def require_columns(row: dict, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in row]
    if missing:
        raise ValueError(f"{label} is missing required fields: {', '.join(missing)}")


def as_bool(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def split_delimited_values(value: str | None) -> list[str]:
    if value is None:
        return []
    parts = [part.strip() for part in value.replace(",", ";").split(";")]
    return [part for part in parts if part and part not in MISSING_VALUE_MARKERS]


def join_values(values: Iterable[str]) -> str:
    cleaned = [value for value in values if value not in MISSING_VALUE_MARKERS]
    return ";".join(sorted(set(cleaned))) if cleaned else NOT_ASSESSED
