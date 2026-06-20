"""Audit the Phase 9 report package for traceability and release readiness."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .schemas import (
    PHASE9_FIGURE_MANIFEST_COLUMNS,
    PHASE9_KEY_FINDINGS_COLUMNS,
    PHASE9_PACKAGE_AUDIT_COLUMNS,
    PHASE9_RELEASE_READINESS_COLUMNS,
)
from .utils import ensure_parent, read_tsv, read_yaml, write_tsv

LOGGER = logging.getLogger(__name__)

FORBIDDEN_REPORT_TERMS = ["activated", "inactivated", "causes", "proves"]
FIGURE_REQUIRED_MARKERS = ["<title", "<desc", "Data source:", "Guardrail:"]


def gitignore_patterns(text: str) -> set[str]:
    return {line.strip() for line in text.splitlines() if line.strip() and not line.strip().startswith("#")}


def release_bundle_paths(config: dict) -> list[str]:
    phase9 = config["phase9_report_generation"]
    package = config.get("phase9_report_package_audit", {})
    if "release_bundle_paths" in package:
        return [Path(path).as_posix() for path in package["release_bundle_paths"]]
    paths = [
        phase9["phase8b_final_evidence"],
        phase9["phase8b_tier_audit"],
        phase9["phase8b_mechanism_summary"],
        phase9["figure_manifest_output"],
        phase9["key_findings_output"],
    ]
    for key in ["package_audit_output", "release_readiness_output"]:
        if key in package:
            paths.append(package[key])
    return [Path(path).as_posix() for path in paths]


def exists_status(path: Path) -> tuple[str, str]:
    if path.exists():
        return "PASS", f"Found {path}"
    return "FAIL", f"Missing {path}"


def semicolon_paths(value: str) -> list[Path]:
    return [Path(part.strip()) for part in value.split(";") if part.strip()]


def add_audit(
    rows: list[dict[str, str]],
    audit_id: str,
    audit_area: str,
    item: str,
    source: str,
    status: str,
    evidence: str,
    required_action: str,
) -> None:
    rows.append(
        {
            "audit_id": audit_id,
            "audit_area": audit_area,
            "item": item,
            "source": source,
            "status": status,
            "evidence": evidence,
            "required_action": required_action,
        }
    )


def add_readiness(
    rows: list[dict[str, str]],
    check_id: str,
    readiness_area: str,
    item: str,
    status: str,
    evidence: str,
    required_action: str,
) -> None:
    rows.append(
        {
            "check_id": check_id,
            "readiness_area": readiness_area,
            "item": item,
            "status": status,
            "evidence": evidence,
            "required_action": required_action,
        }
    )


def audit_report_traceability(config: dict, audit_rows: list[dict[str, str]]) -> None:
    phase9 = config["phase9_report_generation"]
    report = Path(phase9["report_output"])
    status, evidence = exists_status(report)
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-001",
        "report_presence",
        "Final Phase 9 report exists",
        str(report),
        status,
        evidence,
        "Regenerate Phase 9 report if missing.",
    )
    report_text = report.read_text(encoding="utf-8") if report.exists() else ""
    required_inputs = [
        "results/evidence/phase8b_tier_audit.tsv",
        "results/evidence/phase8b_mechanism_summary.tsv",
        "results/evidence/phase8b_final_integrated_evidence.tsv",
    ]
    missing_inputs = [path for path in required_inputs if path not in report_text]
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-002",
        "traceability",
        "Primary evidence inputs are cited in report",
        str(report),
        "PASS" if not missing_inputs else "FAIL",
        "All primary Phase 8b inputs are cited." if not missing_inputs else f"Missing citations: {';'.join(missing_inputs)}",
        "Add missing primary input paths to the report.",
    )
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-003",
        "traceability",
        "Claims register is explicitly referenced",
        str(report),
        "PASS" if "docs/claims_register.md" in report_text else "WARN",
        "Report cites docs/claims_register.md." if "docs/claims_register.md" in report_text else "Report mentions claims register but does not cite docs/claims_register.md explicitly.",
        "Add an explicit claims-register path to the report before public release.",
    )
    found_terms = [term for term in FORBIDDEN_REPORT_TERMS if term in report_text.lower()]
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-004",
        "language_guardrails",
        "Forbidden overclaiming terms are absent from report body",
        str(report),
        "PASS" if not found_terms else "FAIL",
        "No forbidden terms detected." if not found_terms else f"Detected terms: {';'.join(found_terms)}",
        "Rewrite report wording to avoid unsupported pathway-state, gene-state, causal, or proof language.",
    )


def audit_figures(config: dict, audit_rows: list[dict[str, str]]) -> None:
    phase9 = config["phase9_report_generation"]
    manifest_path = Path(phase9["figure_manifest_output"])
    status, evidence = exists_status(manifest_path)
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-005",
        "figure_manifest",
        "Figure manifest exists",
        str(manifest_path),
        status,
        evidence,
        "Regenerate Phase 9 report package if missing.",
    )
    if not manifest_path.exists():
        return
    manifest_rows = read_tsv(manifest_path, PHASE9_FIGURE_MANIFEST_COLUMNS)
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-006",
        "figure_manifest",
        "Expected four Phase 9 figures are registered",
        str(manifest_path),
        "PASS" if len(manifest_rows) == 4 else "FAIL",
        f"Registered figure rows: {len(manifest_rows)}",
        "Regenerate Phase 9 figures and manifest.",
    )
    missing_files: list[str] = []
    missing_markers: list[str] = []
    for row in manifest_rows:
        for key in ["figure_data", "output_svg"]:
            path = Path(row[key])
            if not path.exists():
                missing_files.append(str(path))
        svg_path = Path(row["output_svg"])
        if svg_path.exists():
            svg_text = svg_path.read_text(encoding="utf-8", errors="replace")
            for marker in FIGURE_REQUIRED_MARKERS:
                if marker not in svg_text:
                    missing_markers.append(f"{svg_path}:{marker}")
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-007",
        "figure_files",
        "Figure SVG and figure-data TSV files exist",
        str(manifest_path),
        "PASS" if not missing_files else "FAIL",
        "All manifest-listed figure files exist." if not missing_files else f"Missing files: {';'.join(missing_files)}",
        "Regenerate missing figure files before release.",
    )
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-008",
        "figure_interpretability",
        "SVG figures contain metadata, provenance, and guardrail captions",
        str(manifest_path),
        "PASS" if not missing_markers else "FAIL",
        "All SVGs contain required metadata/provenance markers." if not missing_markers else f"Missing markers: {';'.join(missing_markers)}",
        "Regenerate figures with metadata, source labels, and guardrail captions.",
    )
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-009",
        "figure_stack_policy",
        "Figure stack is suitable for current evidence type",
        str(manifest_path),
        "PASS_WITH_CAVEAT",
        "TSV-backed SVGs are appropriate for categorical evidence audits; standard plotting libraries such as matplotlib/ggplot2 are optional for manuscript styling, not required for this repository stage.",
        "If a manuscript or journal figure package is needed, add an optional matplotlib/R export layer that reads the same figure-data TSVs.",
    )


def audit_key_findings(config: dict, audit_rows: list[dict[str, str]]) -> None:
    phase9 = config["phase9_report_generation"]
    key_findings_path = Path(phase9["key_findings_output"])
    status, evidence = exists_status(key_findings_path)
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-010",
        "key_findings",
        "Key findings table exists",
        str(key_findings_path),
        status,
        evidence,
        "Regenerate Phase 9 key findings if missing.",
    )
    if not key_findings_path.exists():
        return
    key_rows = read_tsv(key_findings_path, PHASE9_KEY_FINDINGS_COLUMNS)
    missing_support: list[str] = []
    for row in key_rows:
        for path in semicolon_paths(row["supporting_files"]):
            if not path.exists():
                missing_support.append(f"{row['finding_id']}:{path}")
    add_audit(
        audit_rows,
        "PHASE9-AUDIT-011",
        "key_findings",
        "Key finding support files exist locally",
        str(key_findings_path),
        "PASS" if not missing_support else "WARN",
        "All supporting files exist locally." if not missing_support else f"Missing or ignored support files: {';'.join(missing_support)}",
        "For public release, either include small support tables, publish them as release artifacts, or document how to regenerate them.",
    )


def audit_public_readiness(config: dict, readiness_rows: list[dict[str, str]]) -> None:
    gitignore = Path(".gitignore")
    text = gitignore.read_text(encoding="utf-8") if gitignore.exists() else ""
    patterns = gitignore_patterns(text)
    add_readiness(
        readiness_rows,
        "PHASE9-READY-001",
        "large_file_hygiene",
        "Raw and intermediate data directories are ignored",
        "PASS" if all(pattern in text for pattern in ["data/raw/", "data/interim/", ".snakemake/"]) else "FAIL",
        ".gitignore excludes data/raw/, data/interim/, and .snakemake/." if gitignore.exists() else ".gitignore is missing.",
        "Keep large raw data and workflow caches out of Git.",
    )
    add_readiness(
        readiness_rows,
        "PHASE9-READY-002",
        "large_file_hygiene",
        "Large sequencing file extensions are ignored",
        "PASS" if all(pattern in text for pattern in ["*.fastq.gz", "*.bam", "*.cram"]) else "FAIL",
        ".gitignore excludes common large sequencing files.",
        "Add missing large-file patterns before public release.",
    )
    add_readiness(
        readiness_rows,
        "PHASE9-READY-003",
        "release_artifacts",
        "Final report and figures are outside ignored generated-report path",
        "PASS" if "reports/generated/" in text and "reports/final/" not in text and "reports/figures/" not in text else "WARN",
        "reports/final/ and reports/figures/ are not globally ignored.",
        "Commit final report/figures only if they are intended release artifacts and small enough for Git.",
    )
    add_readiness(
        readiness_rows,
        "PHASE9-READY-004",
        "release_artifacts",
        "Selected final result TSVs are unignored by exact exception",
        "PASS"
        if (
            "results/*" in patterns
            and "results/" not in patterns
            and "!results/evidence/" in patterns
            and "results/evidence/*" in patterns
            and "!results/reporting/" in patterns
            and "results/reporting/*" in patterns
            and all(f"!{path}" in patterns for path in release_bundle_paths(config))
        )
        else "WARN",
        "Most result outputs remain ignored, while the selected final evidence/reporting TSV release bundle is unignored by exact .gitignore exceptions.",
        "Keep the release bundle narrow. If this check warns, add parent-directory and exact-file .gitignore exceptions for only the selected small TSVs.",
    )
    add_readiness(
        readiness_rows,
        "PHASE9-READY-005",
        "reproducibility",
        "Environment and runbook are present",
        "PASS" if Path("environment.yml").exists() and Path("docs/runbook.md").exists() else "FAIL",
        "environment.yml and docs/runbook.md are present.",
        "Keep direct refresh commands and environment requirements synchronized.",
    )
    add_readiness(
        readiness_rows,
        "PHASE9-READY-006",
        "scientific_traceability",
        "Claims register and study design are present",
        "PASS" if Path("docs/claims_register.md").exists() and Path("docs/study_design.md").exists() else "FAIL",
        "docs/claims_register.md and docs/study_design.md are present.",
        "Do not release final report without claim traceability documents.",
    )
    add_readiness(
        readiness_rows,
        "PHASE9-READY-007",
        "validation",
        "Phase 9 report-generation and package-audit tests are present",
        "PASS"
        if Path("tests/test_phase9_report_generation.py").exists() and Path("tests/test_phase9_package_audit.py").exists()
        else "FAIL",
        "tests/test_phase9_report_generation.py and tests/test_phase9_package_audit.py are present.",
        "Run focused and full tests before public release.",
    )


def write_audit_report(
    audit_rows: list[dict[str, str]],
    readiness_rows: list[dict[str, str]],
    output: Path,
) -> None:
    audit_status_counts = {status: sum(1 for row in audit_rows if row["status"] == status) for status in sorted({row["status"] for row in audit_rows})}
    readiness_status_counts = {
        status: sum(1 for row in readiness_rows if row["status"] == status) for status in sorted({row["status"] for row in readiness_rows})
    }
    action_rows = [row for row in audit_rows if row["status"] != "PASS"]
    readiness_actions = [row for row in readiness_rows if row["status"] != "PASS"]
    lines = [
        "# Phase 9 Report Package Audit",
        "",
        "This audit checks the final Phase 9 report package for traceability, figure interpretability, and public-repository readiness. It does not change evidence tiers or biological interpretation.",
        "",
        "## Figure Stack Decision",
        "",
        "The current TSV-backed SVG approach is suitable for this repository stage because Phase 8b outputs are categorical evidence-audit tables. A standard plotting stack such as matplotlib, seaborn, or ggplot2 can be useful later for manuscript styling, but it is not required for defensible repository figures as long as the figures remain data-backed, accessible, and provenance-labelled.",
        "",
        "## Status Summary",
        "",
        f"- Traceability audit status counts: {audit_status_counts}",
        f"- Release-readiness status counts: {readiness_status_counts}",
        "",
        "## Non-PASS Audit Items",
        "",
    ]
    if action_rows:
        for row in action_rows:
            lines.append(f"- `{row['audit_id']}` {row['status']}: {row['item']}. Action: {row['required_action']}")
    else:
        lines.append("- None.")
    lines.extend(["", "## Non-PASS Release-Readiness Items", ""])
    if readiness_actions:
        for row in readiness_actions:
            lines.append(f"- `{row['check_id']}` {row['status']}: {row['item']}. Action: {row['required_action']}")
    else:
        lines.append("- None.")
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "This audit is about package quality and traceability. It is not evidence that any candidate gene, pathway, repeat context, or expression signal contributes to Greenland shark longevity.",
        ]
    )
    ensure_parent(output)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_audit(config_path: Path, package_audit_output: Path, release_readiness_output: Path, report_output: Path) -> None:
    config = read_yaml(config_path)
    if "phase9_report_generation" not in config:
        raise ValueError("config.yaml is missing phase9_report_generation")
    audit_rows: list[dict[str, str]] = []
    readiness_rows: list[dict[str, str]] = []
    audit_report_traceability(config, audit_rows)
    audit_figures(config, audit_rows)
    audit_key_findings(config, audit_rows)
    audit_public_readiness(config, readiness_rows)
    write_tsv(package_audit_output, audit_rows, PHASE9_PACKAGE_AUDIT_COLUMNS)
    write_tsv(release_readiness_output, readiness_rows, PHASE9_RELEASE_READINESS_COLUMNS)
    write_audit_report(audit_rows, readiness_rows, report_output)
    LOGGER.info("Wrote Phase 9 package audit with %d traceability checks and %d readiness checks", len(audit_rows), len(readiness_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit Phase 9 report package traceability and release readiness.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--package-audit-output", type=Path, required=True)
    parser.add_argument("--release-readiness-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_audit(args.config, args.package_audit_output, args.release_readiness_output, args.report_output)


if __name__ == "__main__":
    main()
