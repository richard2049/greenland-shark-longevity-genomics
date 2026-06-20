"""Run and parse HMMER/Pfam candidate-domain validation."""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .schemas import DOMAIN_CHECK_TARGET_COLUMNS, DOMAIN_INTEGRITY_COLUMNS, HMMER_PFAM_PREFLIGHT_COLUMNS
from .utils import configure_logging, read_tsv, read_yaml, write_tsv


DOMAIN_SUPPORTED = "DOMAIN_SUPPORTED"
PARTIAL_DOMAIN = "PARTIAL_DOMAIN"
NO_EXPECTED_DOMAIN_DETECTED = "NO_EXPECTED_DOMAIN_DETECTED"
NOT_ASSESSED = "NOT_ASSESSED"


@dataclass(frozen=True)
class PfamHit:
    """One parsed HMMER --domtblout domain hit."""

    target_name: str
    target_accession: str
    target_length: int
    query_name: str
    query_accession: str
    query_length: int
    sequence_evalue: float
    sequence_score: float
    domain_number: int
    domain_count: int
    conditional_evalue: float
    independent_evalue: float
    domain_score: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    accuracy: float
    description: str

    @property
    def accession_root(self) -> str:
        return self.target_accession.split(".")[0] if self.target_accession else self.target_name

    @property
    def hmm_coverage(self) -> float:
        if self.target_length <= 0:
            return 0.0
        return max(0, self.hmm_to - self.hmm_from + 1) / self.target_length

    @property
    def query_coverage(self) -> float:
        if self.query_length <= 0:
            return 0.0
        return max(0, self.ali_to - self.ali_from + 1) / self.query_length


def as_float(value: str, label: str) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(f"Cannot parse {label} as float: {value}") from exc


def as_int(value: str, label: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"Cannot parse {label} as integer: {value}") from exc


def parse_domtblout(path: Path) -> list[PfamHit]:
    hits: list[PfamHit] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(maxsplit=22)
            if len(parts) < 22:
                raise ValueError(f"{path}:{line_number} is not a valid HMMER domtblout row")
            description = parts[22] if len(parts) > 22 else ""
            hits.append(
                PfamHit(
                    target_name=parts[0],
                    target_accession=parts[1] if parts[1] != "-" else "",
                    target_length=as_int(parts[2], "target_length"),
                    query_name=parts[3],
                    query_accession=parts[4] if parts[4] != "-" else "",
                    query_length=as_int(parts[5], "query_length"),
                    sequence_evalue=as_float(parts[6], "sequence_evalue"),
                    sequence_score=as_float(parts[7], "sequence_score"),
                    domain_number=as_int(parts[9], "domain_number"),
                    domain_count=as_int(parts[10], "domain_count"),
                    conditional_evalue=as_float(parts[11], "conditional_evalue"),
                    independent_evalue=as_float(parts[12], "independent_evalue"),
                    domain_score=as_float(parts[13], "domain_score"),
                    hmm_from=as_int(parts[15], "hmm_from"),
                    hmm_to=as_int(parts[16], "hmm_to"),
                    ali_from=as_int(parts[17], "ali_from"),
                    ali_to=as_int(parts[18], "ali_to"),
                    env_from=as_int(parts[19], "env_from"),
                    env_to=as_int(parts[20], "env_to"),
                    accuracy=as_float(parts[21], "accuracy"),
                    description=description,
                )
            )
    return hits


def pressed_files_for(hmm_database: Path) -> list[Path]:
    return [Path(str(hmm_database) + suffix) for suffix in [".h3f", ".h3i", ".h3m", ".h3p"]]


def domain_config(config_path: Path) -> dict:
    config = read_yaml(config_path)
    params = config.get("domain_validation", {})
    if not isinstance(params, dict):
        raise ValueError("config/domain_validation must be a mapping")
    return params


def pyhmmer_available() -> tuple[bool, str]:
    try:
        import pyhmmer  # type: ignore
    except ImportError:
        return False, "pyhmmer is not installed in the active Python environment."
    return True, f"pyhmmer {pyhmmer.__version__}"


def build_preflight_rows(config_path: Path, input_fasta: Path) -> list[dict[str, str]]:
    params = domain_config(config_path)
    backend = str(params.get("backend", "hmmscan")).lower()
    executable = str(params.get("hmmer_executable", "hmmscan"))
    hmm_database = Path(str(params.get("pfam_hmm", "NOT_ASSESSED")))
    executable_path = shutil.which(executable)
    has_pyhmmer, pyhmmer_value = pyhmmer_available()

    rows = [
        {
            "check_name": "candidate_fasta",
            "status": "PASS" if input_fasta.exists() else "FAIL",
            "value": str(input_fasta),
            "notes": "Representative candidate FASTA exists." if input_fasta.exists() else "Representative candidate FASTA is missing.",
        },
        {
            "check_name": "domain_backend",
            "status": "PASS" if backend in {"hmmscan", "pyhmmer"} else "FAIL",
            "value": backend,
            "notes": "Supported backend selected." if backend in {"hmmscan", "pyhmmer"} else "Use backend 'hmmscan' or 'pyhmmer'.",
        },
        {
            "check_name": "hmmer_executable",
            "status": "PASS" if backend == "pyhmmer" or executable_path else "FAIL",
            "value": executable_path or executable,
            "notes": "Not required for pyhmmer backend." if backend == "pyhmmer" else ("HMMER hmmscan executable is available on PATH." if executable_path else "Install HMMER or run this step in a container/WSL environment."),
        },
        {
            "check_name": "pyhmmer_package",
            "status": "PASS" if backend == "hmmscan" or has_pyhmmer else "FAIL",
            "value": pyhmmer_value,
            "notes": "Not required for hmmscan backend." if backend == "hmmscan" else ("pyhmmer package is available." if has_pyhmmer else "Install pyhmmer into the active Python environment."),
        },
        {
            "check_name": "pfam_hmm",
            "status": "PASS" if hmm_database.exists() else "FAIL",
            "value": str(hmm_database),
            "notes": "Pfam HMM database exists." if hmm_database.exists() else "Set domain_validation.pfam_hmm to a local Pfam-A.hmm path.",
        },
    ]
    missing_pressed = [path.name for path in pressed_files_for(hmm_database) if not path.exists()]
    rows.append(
        {
            "check_name": "pfam_hmmpress_indexes",
            "status": "PASS" if backend == "pyhmmer" or (hmm_database.exists() and not missing_pressed) else "FAIL",
            "value": ",".join(path.name for path in pressed_files_for(hmm_database)),
            "notes": "Not required for pyhmmer backend." if backend == "pyhmmer" else ("HMMER database indexes exist." if hmm_database.exists() and not missing_pressed else f"Run hmmpress on the Pfam HMM database before hmmscan. Missing: {','.join(missing_pressed) or 'NOT_ASSESSED'}"),
        }
    )
    return rows


def require_preflight_ready(config_path: Path, input_fasta: Path) -> None:
    failures = [row for row in build_preflight_rows(config_path, input_fasta) if row["status"] == "FAIL"]
    if failures:
        messages = [f"{row['check_name']}={row['value']} ({row['notes']})" for row in failures]
        raise SystemExit("HMMER/Pfam preflight failed: " + "; ".join(messages))


def run_hmmscan(config_path: Path, input_fasta: Path, domtblout_output: Path, text_output: Path) -> None:
    params = domain_config(config_path)
    require_preflight_ready(config_path, input_fasta)
    backend = str(params.get("backend", "hmmscan")).lower()
    if backend == "pyhmmer":
        run_pyhmmer_scan(config_path, input_fasta, domtblout_output, text_output)
        return
    executable = str(params.get("hmmer_executable", "hmmscan"))
    hmm_database = Path(str(params["pfam_hmm"]))
    cpu = str(params.get("threads", 1))
    use_cut_ga = bool(params.get("use_cut_ga", True))
    max_i_evalue = str(params.get("max_i_evalue", "1e-5"))

    command = [
        executable,
        "--cpu",
        cpu,
        "--domtblout",
        str(domtblout_output),
        "--noali",
    ]
    if use_cut_ga:
        command.append("--cut_ga")
    else:
        command.extend(["-E", max_i_evalue, "--domE", max_i_evalue])
    command.extend([str(hmm_database), str(input_fasta)])

    domtblout_output.parent.mkdir(parents=True, exist_ok=True)
    text_output.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Running HMMER/Pfam: %s", " ".join(command))
    with text_output.open("w", encoding="utf-8") as stdout_handle:
        subprocess.run(command, check=True, stdout=stdout_handle, stderr=subprocess.STDOUT)


def run_pyhmmer_scan(config_path: Path, input_fasta: Path, domtblout_output: Path, text_output: Path) -> None:
    import pyhmmer
    from pyhmmer import easel, plan7

    params = domain_config(config_path)
    hmm_database = Path(str(params["pfam_hmm"]))
    cpu = int(params.get("threads", 1))
    use_cut_ga = bool(params.get("use_cut_ga", True))
    max_i_evalue = float(params.get("max_i_evalue", 1e-5))

    alphabet = easel.Alphabet.amino()
    with easel.SequenceFile(str(input_fasta), digital=True, alphabet=alphabet) as sequence_file:
        sequences = list(sequence_file)
    if not sequences:
        raise ValueError(f"No protein sequences found in {input_fasta}")

    domtblout_output.parent.mkdir(parents=True, exist_ok=True)
    text_output.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Running pyhmmer/Pfam scan of %s sequences against %s", len(sequences), hmm_database)

    options: dict[str, object] = {"cpus": cpu}
    if use_cut_ga:
        options["bit_cutoffs"] = "gathering"
    else:
        options["E"] = max_i_evalue
        options["domE"] = max_i_evalue

    with plan7.HMMFile(str(hmm_database)) as hmm_file:
        try:
            profiles = hmm_file.optimized_profiles()
            profile_mode = "pressed_optimized_profiles"
        except ValueError:
            profiles = hmm_file
            profile_mode = "plain_hmm_prefetch"
        with domtblout_output.open("wb") as domain_handle, text_output.open("w", encoding="utf-8") as text_handle:
            text_handle.write(f"pyhmmer_version\t{pyhmmer.__version__}\n")
            text_handle.write(f"python\t{sys.executable}\n")
            text_handle.write(f"input_fasta\t{input_fasta}\n")
            text_handle.write(f"pfam_hmm\t{hmm_database}\n")
            text_handle.write(f"sequence_count\t{len(sequences)}\n")
            text_handle.write(f"use_cut_ga\t{use_cut_ga}\n")
            text_handle.write(f"profile_mode\t{profile_mode}\n")
            for index, top_hits in enumerate(pyhmmer.hmmer.hmmscan(sequences, profiles, **options), start=1):
                top_hits.write(domain_handle, format="domains", header=index == 1)
                query_name = top_hits.query.name.decode() if isinstance(top_hits.query.name, bytes) else str(top_hits.query.name)
                text_handle.write(f"query\t{query_name}\treported_hits\t{len(top_hits.reported)}\tincluded_hits\t{len(top_hits.included)}\n")


def hits_by_query(hits: Iterable[PfamHit]) -> dict[str, list[PfamHit]]:
    grouped: dict[str, list[PfamHit]] = {}
    for hit in hits:
        grouped.setdefault(hit.query_name, []).append(hit)
    return grouped


def accepted_hits(hits: list[PfamHit], max_i_evalue: float, min_partial_hmm_coverage: float) -> list[PfamHit]:
    return [
        hit
        for hit in hits
        if hit.independent_evalue <= max_i_evalue and hit.hmm_coverage >= min_partial_hmm_coverage
    ]


def best_hit(hits: list[PfamHit]) -> PfamHit | None:
    if not hits:
        return None
    return sorted(hits, key=lambda hit: (hit.independent_evalue, -hit.domain_score, -hit.hmm_coverage))[0]


def classify_hits(hits: list[PfamHit], max_i_evalue: float, min_partial_hmm_coverage: float, min_full_hmm_coverage: float) -> tuple[str, list[PfamHit], PfamHit | None, str]:
    significant = [hit for hit in hits if hit.independent_evalue <= max_i_evalue]
    accepted = accepted_hits(hits, max_i_evalue, min_partial_hmm_coverage)
    best = best_hit(accepted)
    if best and best.hmm_coverage >= min_full_hmm_coverage:
        return DOMAIN_SUPPORTED, accepted, best, "At least one Pfam domain passed the e-value and full-domain coverage thresholds."
    if significant:
        partial_best = best_hit(significant)
        return PARTIAL_DOMAIN, significant, partial_best, "A significant Pfam hit was detected, but domain HMM coverage did not meet the full-domain threshold."
    return NO_EXPECTED_DOMAIN_DETECTED, [], None, "No Pfam domain passed the configured independent e-value threshold; do not interpret this as gene absence or loss."


def format_float(value: float | None) -> str:
    if value is None:
        return NOT_ASSESSED
    return f"{value:.6g}"


def build_domain_integrity_rows(
    domain_targets_path: Path,
    domtblout_path: Path,
    config_path: Path,
) -> list[dict[str, str]]:
    params = domain_config(config_path)
    max_i_evalue = float(params.get("max_i_evalue", 1e-5))
    min_partial_hmm_coverage = float(params.get("min_partial_hmm_coverage", 0.20))
    min_full_hmm_coverage = float(params.get("min_full_hmm_coverage", 0.70))
    hmm_database = str(params.get("pfam_hmm", NOT_ASSESSED))
    hmmer_command = str(params.get("recommended_command", "hmmscan --cut_ga --domtblout <output> <Pfam-A.hmm> <candidate_fasta>"))
    grouped_hits = hits_by_query(parse_domtblout(domtblout_path))
    target_rows = read_tsv(domain_targets_path, DOMAIN_CHECK_TARGET_COLUMNS)

    output_rows: list[dict[str, str]] = []
    for target in target_rows:
        fasta_id = target["representative_fasta_id"]
        ready = target["domain_check_status"] == "READY_FOR_DOMAIN_SCAN" and fasta_id != NOT_ASSESSED
        if not ready:
            status = NOT_ASSESSED
            accepted: list[PfamHit] = []
            best: PfamHit | None = None
            note = "Candidate sequence was not ready for domain validation."
        else:
            status, accepted, best, note = classify_hits(
                grouped_hits.get(fasta_id, []),
                max_i_evalue=max_i_evalue,
                min_partial_hmm_coverage=min_partial_hmm_coverage,
                min_full_hmm_coverage=min_full_hmm_coverage,
            )

        accepted_accessions = sorted({hit.accession_root for hit in accepted})
        all_hits = grouped_hits.get(fasta_id, [])
        all_hit_labels = [
            f"{hit.accession_root}:{hit.target_name}:iE={format_float(hit.independent_evalue)}:hmm_cov={format_float(hit.hmm_coverage)}"
            for hit in sorted(all_hits, key=lambda item: (item.independent_evalue, -item.domain_score))
        ]
        output_rows.append(
            {
                "mechanism": target["mechanism"],
                "gene_symbol": target["gene_symbol"],
                "resource_id": target["resource_id"],
                "gene_id": target["gene_id"],
                "orthogroup_id": target["orthogroup_id"],
                "representative_protein_id": target["representative_protein_id"],
                "representative_fasta_id": fasta_id,
                "representative_length_aa": target["representative_length_aa"],
                "domain_validation_status": status,
                "best_pfam_accession": best.accession_root if best else NOT_ASSESSED,
                "best_pfam_name": best.target_name if best else NOT_ASSESSED,
                "best_i_evalue": format_float(best.independent_evalue if best else None),
                "best_bitscore": format_float(best.domain_score if best else None),
                "best_hmm_coverage": format_float(best.hmm_coverage if best else None),
                "best_query_coverage": format_float(best.query_coverage if best else None),
                "accepted_domain_count": str(len(accepted)),
                "accepted_domain_accessions": ",".join(accepted_accessions) if accepted_accessions else NOT_ASSESSED,
                "all_domain_hits": ";".join(all_hit_labels) if all_hit_labels else NOT_ASSESSED,
                "hmmer_domtblout": str(domtblout_path),
                "hmmer_database": hmm_database,
                "hmmer_command": hmmer_command,
                "classification_rule": f"DOMAIN_SUPPORTED requires i-E <= {max_i_evalue:g} and HMM coverage >= {min_full_hmm_coverage:g}; PARTIAL_DOMAIN requires significant hit below full coverage; no qualifying hit is not gene absence.",
                "claim_language_guardrail": "Do not infer activation, inactivation, gene absence, functional advantage, or causal longevity relevance from this table alone.",
                "notes": note,
            }
        )
    return output_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="HMMER/Pfam candidate-domain validation.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    preflight_parser = subparsers.add_parser("preflight")
    preflight_parser.add_argument("--config", type=Path, required=True)
    preflight_parser.add_argument("--input-fasta", type=Path, required=True)
    preflight_parser.add_argument("--output", type=Path, required=True)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--config", type=Path, required=True)
    run_parser.add_argument("--input-fasta", type=Path, required=True)
    run_parser.add_argument("--domtblout-output", type=Path, required=True)
    run_parser.add_argument("--text-output", type=Path, required=True)

    parse_parser = subparsers.add_parser("parse")
    parse_parser.add_argument("--config", type=Path, required=True)
    parse_parser.add_argument("--domain-targets", type=Path, required=True)
    parse_parser.add_argument("--domtblout", type=Path, required=True)
    parse_parser.add_argument("--output", type=Path, required=True)

    args = parser.parse_args()
    configure_logging()
    if args.command == "preflight":
        write_tsv(args.output, build_preflight_rows(args.config, args.input_fasta), HMMER_PFAM_PREFLIGHT_COLUMNS)
    elif args.command == "run":
        run_hmmscan(args.config, args.input_fasta, args.domtblout_output, args.text_output)
    elif args.command == "parse":
        rows = build_domain_integrity_rows(args.domain_targets, args.domtblout, args.config)
        write_tsv(args.output, rows, DOMAIN_INTEGRITY_COLUMNS)
    else:
        raise ValueError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
