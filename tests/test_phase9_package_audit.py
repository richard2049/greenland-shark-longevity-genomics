from pathlib import Path

import yaml

from greenland_shark_longevity.phase9_package_audit import run_audit
from greenland_shark_longevity.schemas import PHASE9_PACKAGE_AUDIT_COLUMNS, PHASE9_RELEASE_READINESS_COLUMNS
from greenland_shark_longevity.utils import read_tsv, write_tsv


def write_svg(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                '<svg xmlns="http://www.w3.org/2000/svg">',
                "<title>Figure</title>",
                "<desc>Test figure.</desc>",
                '<text x="1" y="1">Data source: test.tsv</text>',
                '<text x="1" y="2">Guardrail: reporting only</text>',
                "</svg>",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_phase9_package_audit_reports_traceability_and_release_bundle_status(tmp_path: Path) -> None:
    report = tmp_path / "reports" / "final" / "phase9.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "\n".join(
            [
                "# Report",
                "results/evidence/phase8b_tier_audit.tsv",
                "results/evidence/phase8b_mechanism_summary.tsv",
                "results/evidence/phase8b_final_integrated_evidence.tsv",
                "docs/claims_register.md",
                "No pathway-state claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "results" / "reporting" / "phase9_figure_manifest.tsv"
    figure_rows = []
    for index in range(4):
        data = tmp_path / "reports" / "figures" / "data" / f"figure_{index}.tsv"
        svg = tmp_path / "reports" / "figures" / f"figure_{index}.svg"
        write_tsv(data, [{"category": "a", "count": "1"}], ["category", "count"])
        write_svg(svg)
        figure_rows.append(
            {
                "figure_id": f"FIG9-{index + 1:03d}",
                "title": "Figure",
                "figure_type": "horizontal_bar_svg",
                "source_table": "test.tsv",
                "figure_data": str(data),
                "output_svg": str(svg),
                "status": "GENERATED",
                "interpretation_guardrail": "Reporting only.",
            }
        )
    write_tsv(
        manifest,
        figure_rows,
        ["figure_id", "title", "figure_type", "source_table", "figure_data", "output_svg", "status", "interpretation_guardrail"],
    )
    support = tmp_path / "support.tsv"
    write_tsv(support, [{"x": "1"}], ["x"])
    key_findings = tmp_path / "results" / "reporting" / "phase9_key_findings.tsv"
    write_tsv(
        key_findings,
        [
            {
                "finding_id": "KF",
                "question": "Q",
                "answer": "A",
                "supporting_files": str(support),
                "interpretation_guardrail": "Reporting only.",
            }
        ],
        ["finding_id", "question", "answer", "supporting_files", "interpretation_guardrail"],
    )
    config = {
        "phase9_report_generation": {
            "phase8b_final_evidence": "results/evidence/phase8b_final_integrated_evidence.tsv",
            "phase8b_tier_audit": "results/evidence/phase8b_tier_audit.tsv",
            "phase8b_mechanism_summary": "results/evidence/phase8b_mechanism_summary.tsv",
            "report_output": str(report),
            "figure_manifest_output": str(manifest),
            "key_findings_output": str(key_findings),
        },
        "phase9_report_package_audit": {
            "release_bundle_paths": [
                "results/evidence/phase8b_final_integrated_evidence.tsv",
                "results/evidence/phase8b_tier_audit.tsv",
                "results/evidence/phase8b_mechanism_summary.tsv",
                "results/reporting/phase9_figure_manifest.tsv",
                "results/reporting/phase9_key_findings.tsv",
                "results/reporting/phase9_report_package_audit.tsv",
                "results/reporting/phase9_public_repository_readiness.tsv",
            ],
            "package_audit_output": "results/reporting/phase9_report_package_audit.tsv",
            "release_readiness_output": "results/reporting/phase9_public_repository_readiness.tsv",
        }
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    package_audit = tmp_path / "results" / "reporting" / "phase9_report_package_audit.tsv"
    readiness = tmp_path / "results" / "reporting" / "phase9_public_repository_readiness.tsv"
    audit_report = tmp_path / "reports" / "final" / "phase9_report_package_audit.md"

    run_audit(config_path, package_audit, readiness, audit_report)

    audit_rows = read_tsv(package_audit, PHASE9_PACKAGE_AUDIT_COLUMNS)
    readiness_rows = read_tsv(readiness, PHASE9_RELEASE_READINESS_COLUMNS)
    assert any(row["audit_id"] == "PHASE9-AUDIT-009" and row["status"] == "PASS_WITH_CAVEAT" for row in audit_rows)
    assert any(row["audit_id"] == "PHASE9-AUDIT-008" and row["status"] == "PASS" for row in audit_rows)
    assert any(row["check_id"] == "PHASE9-READY-004" and row["status"] == "PASS" for row in readiness_rows)
    assert "Figure Stack Decision" in audit_report.read_text(encoding="utf-8")
