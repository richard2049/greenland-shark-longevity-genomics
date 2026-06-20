"""OrthoFinder preflight checks and result parsing."""

from __future__ import annotations

import argparse
import csv
import platform
import shutil
import subprocess
from pathlib import Path

from .schemas import (
    ORTHOFINDER_INPUT_MANIFEST_COLUMNS,
    ORTHOFINDER_PREFLIGHT_COLUMNS,
    ORTHOFINDER_SPECIES_SUMMARY_COLUMNS,
    ORTHOGROUP_GENE_COUNT_LONG_COLUMNS,
)
from .utils import configure_logging, read_tsv, read_yaml, write_tsv


def ready_orthofinder_inputs(manifest_path: Path) -> list[dict[str, str]]:
    rows = read_tsv(manifest_path, ORTHOFINDER_INPUT_MANIFEST_COLUMNS)
    return [row for row in rows if row["orthofinder_ready"] == "True"]


def missing_orthofinder_message(executable: str) -> str:
    if platform.system() == "Windows":
        return (
            f"{executable!r} is not available in the current native Windows environment. "
            "The Bioconda OrthoFinder package depends on external search tools such as BLAST/DIAMOND; "
            "this repo expects the orthology step to run in WSL/Linux or another environment where "
            "OrthoFinder is on PATH."
        )
    return (
        f"{executable!r} is not available on PATH. Install OrthoFinder in the active environment "
        "or set orthofinder.executable in config/config.yaml to the executable path."
    )


def build_preflight_rows(config_path: Path, manifest_path: Path) -> list[dict[str, str]]:
    config = read_yaml(config_path)
    orthofinder_config = config.get("orthofinder", {})
    executable = str(orthofinder_config.get("executable", "orthofinder"))
    minimum_species = int(orthofinder_config.get("minimum_species", 2))
    ready_rows = ready_orthofinder_inputs(manifest_path)
    species_ids = sorted({row["species_id"] for row in ready_rows})
    executable_path = shutil.which(executable)

    rows = [
        {
            "check_name": "orthofinder_binary",
            "status": "PASS" if executable_path else "FAIL",
            "value": executable_path or executable,
            "notes": "OrthoFinder executable found on PATH."
            if executable_path
            else missing_orthofinder_message(executable),
        },
        {
            "check_name": "ready_species_count",
            "status": "PASS" if len(species_ids) >= minimum_species else "FAIL",
            "value": str(len(species_ids)),
            "notes": f"Minimum configured species count is {minimum_species}.",
        },
    ]

    for row in ready_rows:
        path = Path(row["orthofinder_input_path"])
        if not path.exists():
            status = "FAIL"
            notes = "Staged OrthoFinder FASTA is missing."
        elif path.stat().st_size == 0:
            status = "FAIL"
            notes = "Staged OrthoFinder FASTA is empty."
        else:
            status = "PASS"
            notes = "Staged OrthoFinder FASTA exists."
        rows.append(
            {
                "check_name": f"input_fasta:{row['species_id']}:{row['resource_id']}",
                "status": status,
                "value": str(path),
                "notes": notes,
            }
        )
    return rows


def require_preflight(config_path: Path, manifest_path: Path) -> None:
    failures = [row for row in build_preflight_rows(config_path, manifest_path) if row["status"] != "PASS"]
    if failures:
        details = "; ".join(f"{row['check_name']}={row['notes']}" for row in failures)
        raise SystemExit(f"OrthoFinder preflight failed: {details}")


def find_latest_orthofinder_results(output_dir: Path) -> Path:
    direct = output_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    if direct.exists():
        return output_dir
    candidates = [
        path.parent.parent
        for path in output_dir.glob("**/Orthogroups/Orthogroups.GeneCount.tsv")
        if path.is_file()
    ]
    if not candidates:
        raise FileNotFoundError(f"No Orthogroups.GeneCount.tsv found under {output_dir}")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def _column_to_manifest_row(manifest_rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    for row in manifest_rows:
        if row["orthofinder_ready"] != "True":
            continue
        path = Path(row["orthofinder_input_path"])
        keys = {path.name, path.stem, path.name.replace(".faa", ""), path.stem.replace(".faa", "")}
        for key in keys:
            mapping[key] = row
    return mapping


def parse_gene_count_matrix(
    results_dir: Path,
    manifest_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], Path]:
    result_root = find_latest_orthofinder_results(results_dir)
    gene_count_path = result_root / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    manifest_rows = read_tsv(manifest_path, ORTHOFINDER_INPUT_MANIFEST_COLUMNS)
    column_map = _column_to_manifest_row(manifest_rows)
    long_rows: list[dict[str, str]] = []
    summary: dict[str, dict[str, int | str]] = {}

    with gene_count_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "Orthogroup" not in reader.fieldnames:
            raise ValueError(f"{gene_count_path} is missing an Orthogroup column")
        species_columns = [column for column in reader.fieldnames if column not in {"Orthogroup", "Total"}]
        for source_column in species_columns:
            if source_column not in column_map:
                raise ValueError(
                    f"OrthoFinder species column {source_column!r} does not match staged input manifest"
                )
            manifest_row = column_map[source_column]
            summary[source_column] = {
                "species_id": manifest_row["species_id"],
                "resource_id": manifest_row["resource_id"],
                "orthogroups_with_genes": 0,
                "total_gene_copies_in_orthogroups": 0,
            }
        for row in reader:
            orthogroup_id = row["Orthogroup"]
            for source_column in species_columns:
                count = int(row[source_column])
                manifest_row = column_map[source_column]
                if count > 0:
                    summary[source_column]["orthogroups_with_genes"] = int(
                        summary[source_column]["orthogroups_with_genes"]
                    ) + 1
                    summary[source_column]["total_gene_copies_in_orthogroups"] = int(
                        summary[source_column]["total_gene_copies_in_orthogroups"]
                    ) + count
                long_rows.append(
                    {
                        "orthogroup_id": orthogroup_id,
                        "species_id": manifest_row["species_id"],
                        "resource_id": manifest_row["resource_id"],
                        "orthofinder_species_column": source_column,
                        "copy_count": str(count),
                        "orthofinder_results_dir": str(result_root),
                    }
                )

    summary_rows = [
        {
            "species_id": str(values["species_id"]),
            "resource_id": str(values["resource_id"]),
            "orthofinder_species_column": source_column,
            "orthogroups_with_genes": str(values["orthogroups_with_genes"]),
            "total_gene_copies_in_orthogroups": str(values["total_gene_copies_in_orthogroups"]),
            "orthofinder_results_dir": str(result_root),
        }
        for source_column, values in sorted(summary.items())
    ]
    return long_rows, summary_rows, result_root


def run_orthofinder(config_path: Path, manifest_path: Path, marker_output: Path | None = None) -> None:
    require_preflight(config_path, manifest_path)
    config = read_yaml(config_path)
    orthofinder_config = config.get("orthofinder", {})
    executable = str(orthofinder_config.get("executable", "orthofinder"))
    input_dir = str(orthofinder_config.get("input_dir", "data/interim/orthofinder_input"))
    output_dir = str(orthofinder_config.get("output_dir", "results/orthofinder"))
    threads = str(orthofinder_config.get("threads", 4))
    command = [executable, "-f", input_dir, "-o", output_dir, "-t", threads, "-a", threads]
    subprocess.run(command, check=True)
    if marker_output is not None:
        marker_output.parent.mkdir(parents=True, exist_ok=True)
        marker_output.write_text("completed\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run OrthoFinder preflight checks or parse results.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    preflight_parser = subparsers.add_parser("preflight")
    preflight_parser.add_argument("--config", type=Path, required=True)
    preflight_parser.add_argument("--manifest", type=Path, required=True)
    preflight_parser.add_argument("--output", type=Path, required=True)

    require_parser = subparsers.add_parser("require")
    require_parser.add_argument("--config", type=Path, required=True)
    require_parser.add_argument("--manifest", type=Path, required=True)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--config", type=Path, required=True)
    run_parser.add_argument("--manifest", type=Path, required=True)
    run_parser.add_argument("--marker-output", type=Path)

    parse_parser = subparsers.add_parser("parse")
    parse_parser.add_argument("--results-dir", type=Path, required=True)
    parse_parser.add_argument("--manifest", type=Path, required=True)
    parse_parser.add_argument("--gene-count-output", type=Path, required=True)
    parse_parser.add_argument("--species-summary-output", type=Path, required=True)

    args = parser.parse_args()
    configure_logging()

    if args.command == "preflight":
        write_tsv(args.output, build_preflight_rows(args.config, args.manifest), ORTHOFINDER_PREFLIGHT_COLUMNS)
    elif args.command == "require":
        require_preflight(args.config, args.manifest)
    elif args.command == "run":
        run_orthofinder(args.config, args.manifest, args.marker_output)
    elif args.command == "parse":
        long_rows, summary_rows, _result_root = parse_gene_count_matrix(args.results_dir, args.manifest)
        write_tsv(args.gene_count_output, long_rows, ORTHOGROUP_GENE_COUNT_LONG_COLUMNS)
        write_tsv(args.species_summary_output, summary_rows, ORTHOFINDER_SPECIES_SUMMARY_COLUMNS)


if __name__ == "__main__":
    main()
