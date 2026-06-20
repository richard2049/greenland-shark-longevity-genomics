"""Phase 9 report and figure generation from final evidence audit tables."""

from __future__ import annotations

import argparse
import html
import logging
from collections import Counter
from pathlib import Path
from textwrap import shorten

from .schemas import (
    PHASE8B_MECHANISM_SUMMARY_COLUMNS,
    PHASE8B_TIER_AUDIT_COLUMNS,
    PHASE9_FIGURE_MANIFEST_COLUMNS,
    PHASE9_KEY_FINDINGS_COLUMNS,
)
from .utils import ensure_parent, read_tsv, read_yaml, write_tsv

LOGGER = logging.getLogger(__name__)

TIER_ORDER = ["Tier 1", "Tier 2", "Tier 3", "Artifact/uncertain"]

# Okabe-Ito inspired colors keep evidence classes distinguishable for common
# color-vision deficiencies while remaining readable in print and slides.
TIER_COLORS = {
    "Tier 1": "#009E73",
    "Tier 2": "#0072B2",
    "Tier 3": "#E69F00",
    "Artifact/uncertain": "#D55E00",
}
CLASS_COLORS = {
    "ROBUST_COMPUTATIONAL_EVIDENCE_READY_FOR_REPORTING": "#009E73",
    "PLAUSIBLE_LEAD_REQUIRES_VALIDATION": "#0072B2",
    "PLAUSIBLE_LEAD_ARTIFACT_PRONE": "#E69F00",
    "ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY": "#D55E00",
    "EXPLORATORY_SIGNAL": "#999999",
    "EXPLORATORY_SIGNAL_WITH_RETINA_CONTEXT": "#56B4E9",
}
CLASS_LABELS = {
    "ROBUST_COMPUTATIONAL_EVIDENCE_READY_FOR_REPORTING": "Robust computational evidence",
    "PLAUSIBLE_LEAD_REQUIRES_VALIDATION": "Plausible lead, requires validation",
    "PLAUSIBLE_LEAD_ARTIFACT_PRONE": "Plausible lead, artifact-prone",
    "ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY": "Artifact/uncertain, do not claim",
    "EXPLORATORY_SIGNAL": "Exploratory signal",
    "EXPLORATORY_SIGNAL_WITH_RETINA_CONTEXT": "Exploratory retina-context signal",
}
CLASS_ORDER = [
    "ROBUST_COMPUTATIONAL_EVIDENCE_READY_FOR_REPORTING",
    "PLAUSIBLE_LEAD_REQUIRES_VALIDATION",
    "PLAUSIBLE_LEAD_ARTIFACT_PRONE",
    "ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY",
    "EXPLORATORY_SIGNAL_WITH_RETINA_CONTEXT",
    "EXPLORATORY_SIGNAL",
]
ARTIFACT_LABELS = {
    "artifact_risk_high": "High artifact risk",
    "artifact_risk_moderate": "Moderate artifact risk",
    "repeat_context_caveat": "Repeat-context caveat",
    "locus_review_required": "Locus review required",
    "retina_expression_supported": "Retina expression support",
}
ARTIFACT_COLORS = {
    "artifact_risk_high": "#D55E00",
    "artifact_risk_moderate": "#E69F00",
    "repeat_context_caveat": "#CC79A7",
    "locus_review_required": "#999999",
    "retina_expression_supported": "#56B4E9",
}

TIER_SUMMARY_COLUMNS = ["evidence_tier", "count", "fraction"]
CLASS_SUMMARY_COLUMNS = ["reporting_class", "count", "fraction"]
MECHANISM_MATRIX_COLUMNS = [
    "mechanism",
    "candidate_count",
    "tier1_count",
    "tier2_count",
    "tier3_count",
    "artifact_uncertain_count",
    "plausible_lead_count",
    "artifact_prone_count",
    "expression_supported_count",
    "repeat_context_caveat_count",
    "locus_review_required_count",
]
ARTIFACT_CONTEXT_COLUMNS = ["category", "count", "source", "interpretation"]


def to_int(value: str | int | None) -> int:
    if value in {None, ""}:
        return 0
    return int(float(str(value)))


def fraction(count: int, total: int) -> str:
    if total == 0:
        return "0.000"
    return f"{count / total:.3f}"


def svg_text(
    x: int,
    y: int,
    text: str,
    size: int = 12,
    weight: str = "400",
    fill: str = "#222222",
    anchor: str = "start",
) -> str:
    return (
        f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">'
        f"{html.escape(text)}</text>"
    )


def svg_multiline_text(
    x: int,
    y: int,
    text: str,
    width: int,
    line_height: int = 14,
    size: int = 12,
    weight: str = "400",
    fill: str = "#222222",
) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current: list[str] = []
    max_chars = max(10, int(width / max(size * 0.58, 1)))
    for word in words:
        candidate = " ".join(current + [word])
        if len(candidate) > max_chars and current:
            lines.append(" ".join(current))
            current = [word]
        else:
            current.append(word)
    if current:
        lines.append(" ".join(current))
    return [
        svg_text(x, y + index * line_height, line, size=size, weight=weight, fill=fill)
        for index, line in enumerate(lines[:3])
    ]


def svg_header(width: int, height: int, title: str, description: str) -> list[str]:
    return [
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">'
        ),
        f"<title id=\"title\">{html.escape(title)}</title>",
        f"<desc id=\"desc\">{html.escape(description)}</desc>",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]


def count_ticks(max_value: int, tick_count: int = 4) -> list[int]:
    if max_value <= 0:
        return [0]
    raw_step = max(1, max_value / tick_count)
    for step in [1, 2, 5, 10, 20, 50, 100]:
        if raw_step <= step:
            chosen = step
            break
    else:
        chosen = int(raw_step)
    ticks = [tick for tick in range(0, max_value + 1, chosen) if tick <= max_value]
    if not ticks or ticks[-1] < max_value:
        ticks.append(max_value)
    return sorted(set(ticks))


def draw_x_axis(left: int, top: int, plot_width: int, plot_height: int, max_value: int) -> list[str]:
    elements = [
        f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#444444" stroke-width="1"/>',
    ]
    for tick in count_ticks(max_value):
        x = left + int((tick / max(max_value, 1)) * plot_width)
        elements.extend(
            [
                f'<line x1="{x}" y1="{top}" x2="{x}" y2="{top + plot_height}" stroke="#e3e3e3" stroke-width="1"/>',
                f'<line x1="{x}" y1="{top + plot_height}" x2="{x}" y2="{top + plot_height + 5}" stroke="#444444" stroke-width="1"/>',
                svg_text(x, top + plot_height + 20, str(tick), size=11, fill="#444444", anchor="middle"),
            ]
        )
    elements.append(svg_text(left + plot_width // 2, top + plot_height + 42, "Candidate count", size=12, fill="#444444", anchor="middle"))
    return elements


def legend_item(x: int, y: int, label: str, color: str) -> list[str]:
    return [
        f'<rect x="{x}" y="{y - 10}" width="16" height="12" rx="1" ry="1" fill="{color}"/>',
        svg_text(x + 22, y, label, size=12, fill="#333333"),
    ]


def write_horizontal_bar_svg(
    rows: list[dict[str, str]],
    label_key: str,
    value_key: str,
    output: Path,
    title: str,
    subtitle: str,
    color_map: dict[str, str] | None = None,
    label_map: dict[str, str] | None = None,
    source_label: str = "Phase 8b audit table",
    guardrail: str = "Counts summarize reporting readiness only.",
    width: int = 1050,
) -> None:
    if not rows:
        raise ValueError(f"No rows available for {title}")
    left = 370
    right = 110
    top = 116
    row_height = 44
    bar_height = 22
    max_value = max(to_int(row[value_key]) for row in rows) or 1
    plot_width = width - left - right
    plot_height = len(rows) * row_height
    height = top + plot_height + 96
    elements = [
        *svg_header(width, height, title, subtitle),
        svg_text(24, 34, title, size=20, weight="700"),
        svg_text(24, 58, subtitle, size=12, fill="#555555"),
        svg_text(24, 84, "Source values are written to the matching reports/figures/data TSV.", size=11, fill="#666666"),
        *draw_x_axis(left, top - 8, plot_width, plot_height + 8, max_value),
    ]
    for index, row in enumerate(rows):
        label = row[label_key]
        display_label = label_map.get(label, label) if label_map else label
        value = to_int(row[value_key])
        percent = row.get("fraction", "")
        y = top + index * row_height
        bar_width = int((value / max_value) * plot_width)
        color = color_map.get(label, "#4E79A7") if color_map else "#4E79A7"
        elements.extend(
            [
                f'<line x1="{left}" y1="{y + bar_height // 2}" x2="{left - 8}" y2="{y + bar_height // 2}" stroke="#444444" stroke-width="1"/>',
                *svg_multiline_text(24, y + 13, display_label, width=320, line_height=13, size=12),
                f'<rect x="{left}" y="{y}" width="{plot_width}" height="{bar_height}" fill="#eeeeee"/>',
                f'<rect x="{left}" y="{y}" width="{bar_width}" height="{bar_height}" fill="{color}"/>',
                svg_text(left + plot_width + 14, y + 16, str(value), size=12, weight="700"),
            ]
        )
        if percent:
            elements.append(svg_text(left + plot_width + 48, y + 16, f"({float(percent) * 100:.1f}%)", size=11, fill="#555555"))
    elements.extend(
        [
            svg_text(24, height - 34, f"Data source: {source_label}", size=11, fill="#555555"),
            svg_text(24, height - 16, f"Guardrail: {guardrail}", size=11, fill="#555555"),
            "</svg>",
        ]
    )
    ensure_parent(output)
    output.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_mechanism_matrix_svg(rows: list[dict[str, str]], output: Path) -> None:
    if not rows:
        raise ValueError("No mechanism rows available for Phase 9 matrix")
    left = 380
    plot_width = 520
    context_left = left + plot_width + 96
    top = 138
    row_height = 46
    width = 1320
    height = top + len(rows) * row_height + 124
    plot_height = len(rows) * row_height
    max_total = max(to_int(row["candidate_count"]) for row in rows) or 1
    elements = [
        *svg_header(
            width,
            height,
            "Phase 9 mechanism evidence matrix",
            "Stacked Phase 8b evidence-tier counts with separate context counters for expression, repeats, and locus review.",
        ),
        svg_text(24, 34, "Phase 9 mechanism evidence matrix", size=20, weight="700"),
        svg_text(
            24,
            58,
            "Stacked counts from Phase 8b mechanism summary; artifact/context counts are not mechanism claims.",
            size=12,
            fill="#555555",
        ),
        svg_text(24, 84, "Source values are written to reports/figures/data/phase9_mechanism_evidence_matrix.tsv.", size=11, fill="#666666"),
        *legend_item(left, 106, "Tier 1", TIER_COLORS["Tier 1"]),
        *legend_item(left + 84, 106, "Tier 2", TIER_COLORS["Tier 2"]),
        *legend_item(left + 168, 106, "Tier 3", TIER_COLORS["Tier 3"]),
        *legend_item(left + 252, 106, "Artifact/uncertain", TIER_COLORS["Artifact/uncertain"]),
        svg_text(context_left, 106, "Context counters", size=12, weight="700", fill="#333333"),
        svg_text(context_left, 128, "expr", size=11, fill="#333333"),
        svg_text(context_left + 72, 128, "repeat", size=11, fill="#333333"),
        svg_text(context_left + 154, 128, "locus", size=11, fill="#333333"),
        *draw_x_axis(left, top - 10, plot_width, plot_height + 10, max_total),
    ]
    for index, row in enumerate(rows):
        y = top + index * row_height
        tier_values = [
            ("Tier 1", to_int(row["tier1_count"])),
            ("Tier 2", to_int(row["tier2_count"])),
            ("Tier 3", to_int(row["tier3_count"])),
            ("Artifact/uncertain", to_int(row["artifact_uncertain_count"])),
        ]
        total = to_int(row["candidate_count"])
        expression = to_int(row["expression_supported_count"])
        repeat = to_int(row["repeat_context_caveat_count"])
        locus = to_int(row["locus_review_required_count"])
        current_x = left
        elements.extend(
            [
                *svg_multiline_text(24, y + 14, row["mechanism"], width=330, line_height=13, size=12),
                f'<rect x="{left}" y="{y}" width="{plot_width}" height="24" fill="#eeeeee"/>',
            ]
        )
        for tier, count in tier_values:
            segment_width = int((count / max_total) * plot_width)
            if count > 0:
                elements.append(
                    f'<rect x="{current_x}" y="{y}" width="{segment_width}" height="24" fill="{TIER_COLORS[tier]}"/>'
                )
            current_x += segment_width
        elements.extend(
            [
                svg_text(left + plot_width + 12, y + 17, f"n={total}", size=12, weight="700"),
                f'<circle cx="{context_left + 12}" cy="{y + 12}" r="8" fill="#56B4E9"/>',
                svg_text(context_left + 12, y + 16, str(expression), size=10, weight="700", fill="#ffffff", anchor="middle"),
                f'<circle cx="{context_left + 92}" cy="{y + 12}" r="8" fill="#CC79A7"/>',
                svg_text(context_left + 92, y + 16, str(repeat), size=10, weight="700", fill="#ffffff", anchor="middle"),
                f'<circle cx="{context_left + 170}" cy="{y + 12}" r="8" fill="#999999"/>',
                svg_text(context_left + 170, y + 16, str(locus), size=10, weight="700", fill="#ffffff", anchor="middle"),
            ]
        )
    elements.extend(
        [
            svg_text(
                24,
                height - 42,
                "Data source: results/evidence/phase8b_mechanism_summary.tsv",
                size=11,
                fill="#555555",
            ),
            svg_text(
                24,
                height - 22,
                "Guardrail: mechanism-level bars summarize reporting readiness only; row-level validation controls biological wording.",
                size=11,
                fill="#555555",
            ),
            "</svg>",
        ]
    )
    ensure_parent(output)
    output.write_text("\n".join(elements) + "\n", encoding="utf-8")


def build_tier_summary(audit_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    counts = Counter(row["phase8b_final_evidence_tier"] for row in audit_rows)
    total = len(audit_rows)
    return [
        {"evidence_tier": tier, "count": str(counts.get(tier, 0)), "fraction": fraction(counts.get(tier, 0), total)}
        for tier in TIER_ORDER
    ]


def build_class_summary(audit_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    counts = Counter(row["phase8b_reporting_class"] for row in audit_rows)
    total = len(audit_rows)
    labels = [label for label in CLASS_ORDER if label in counts]
    labels.extend(sorted(label for label in counts if label not in CLASS_ORDER))
    return [{"reporting_class": label, "count": str(counts[label]), "fraction": fraction(counts[label], total)} for label in labels]


def build_mechanism_matrix(summary_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in summary_rows:
        rows.append({column: row.get(column, "0") for column in MECHANISM_MATRIX_COLUMNS})
    rows.sort(key=lambda item: (to_int(item["artifact_uncertain_count"]), to_int(item["artifact_prone_count"]), item["mechanism"]), reverse=True)
    return rows


def build_artifact_context_rows(audit_rows: list[dict[str, str]], summary_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    risk_counts = Counter(row["artifact_risk_level"] for row in audit_rows)
    caveat_rows = [
        {
            "category": "artifact_risk_high",
            "count": str(risk_counts.get("high", 0)),
            "source": "phase8b_tier_audit.tsv:artifact_risk_level",
            "interpretation": "High artifact risk rows require row-level review before biological wording.",
        },
        {
            "category": "artifact_risk_moderate",
            "count": str(risk_counts.get("moderate", 0)),
            "source": "phase8b_tier_audit.tsv:artifact_risk_level",
            "interpretation": "Moderate artifact risk rows remain plausible only when validation caveats are retained.",
        },
        {
            "category": "repeat_context_caveat",
            "count": str(sum(to_int(row["repeat_context_caveat_count"]) for row in summary_rows)),
            "source": "phase8b_mechanism_summary.tsv:repeat_context_caveat_count",
            "interpretation": "Repeat context is artifact/context evidence only.",
        },
        {
            "category": "locus_review_required",
            "count": str(sum(to_int(row["locus_review_required_count"]) for row in summary_rows)),
            "source": "phase8b_mechanism_summary.tsv:locus_review_required_count",
            "interpretation": "Locus review is required before duplication or paralog-specific wording.",
        },
        {
            "category": "retina_expression_supported",
            "count": str(sum(to_int(row["expression_supported_count"]) for row in summary_rows)),
            "source": "phase8b_mechanism_summary.tsv:expression_supported_count",
            "interpretation": "Expression support is retina-specific and not a pathway-state result.",
        },
    ]
    return caveat_rows


def build_key_findings(audit_rows: list[dict[str, str]], summary_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    tiers = Counter(row["phase8b_final_evidence_tier"] for row in audit_rows)
    classes = Counter(row["phase8b_reporting_class"] for row in audit_rows)
    artifact_prone_genes = sorted(
        row["gene_symbol"]
        for row in audit_rows
        if row["phase8b_reporting_class"] == "PLAUSIBLE_LEAD_ARTIFACT_PRONE"
    )
    uncertain_genes = sorted(
        row["gene_symbol"]
        for row in audit_rows
        if row["phase8b_reporting_class"] == "ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY"
    )
    expression_supported = sum(to_int(row["expression_supported_count"]) for row in summary_rows)
    return [
        {
            "finding_id": "PHASE9-KF-001",
            "question": "Does any candidate meet current Tier 1 criteria?",
            "answer": f"No current candidate reaches Tier 1 in Phase 8b ({tiers.get('Tier 1', 0)} rows).",
            "supporting_files": "results/evidence/phase8b_tier_audit.tsv;results/evidence/phase8b_mechanism_summary.tsv",
            "interpretation_guardrail": "Report this as current evidence status, not as a statement about biology.",
        },
        {
            "finding_id": "PHASE9-KF-002",
            "question": "What is plausible but incomplete?",
            "answer": f"{tiers.get('Tier 2', 0)} rows remain Tier 2; {classes.get('PLAUSIBLE_LEAD_REQUIRES_VALIDATION', 0)} are plausible leads requiring validation.",
            "supporting_files": "results/evidence/phase8b_tier_audit.tsv",
            "interpretation_guardrail": "Tier 2 supports conservative candidate language only.",
        },
        {
            "finding_id": "PHASE9-KF-003",
            "question": "Which plausible leads are artifact-prone?",
            "answer": ", ".join(artifact_prone_genes) if artifact_prone_genes else "No Phase 8b plausible lead is marked artifact-prone.",
            "supporting_files": "results/evidence/phase8b_tier_audit.tsv;results/repeats/phase5c_gene_repeat_qc_summary.tsv;results/rescue/phase4e_gene_hardened_summary.tsv",
            "interpretation_guardrail": "Artifact-prone leads require locus, repeat, paralog, and cross-resource validation before stronger wording.",
        },
        {
            "finding_id": "PHASE9-KF-004",
            "question": "What should not be claimed biologically?",
            "answer": f"{classes.get('ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY', 0)} rows remain artifact/uncertain: {', '.join(uncertain_genes)}.",
            "supporting_files": "results/evidence/phase8b_tier_audit.tsv",
            "interpretation_guardrail": "Do not turn artifact/uncertain rows into mechanism statements.",
        },
        {
            "finding_id": "PHASE9-KF-005",
            "question": "How should expression support be interpreted?",
            "answer": f"{expression_supported} candidate rows carry cautious retina-specific expression support.",
            "supporting_files": "results/evidence/phase8b_mechanism_summary.tsv;results/rnaseq/phase7e_candidate_expression_hardened.tsv",
            "interpretation_guardrail": "Retina expression support is tissue-specific and not differential expression.",
        },
    ]


def relative_for_markdown(target: Path, from_file: Path) -> str:
    try:
        return target.resolve().relative_to(from_file.parent.resolve()).as_posix()
    except ValueError:
        import os

        return Path(os.path.relpath(target.resolve(), from_file.parent.resolve())).as_posix()


def path_text(path: Path) -> str:
    return path.as_posix()


def markdown_table(rows: list[dict[str, str]], columns: list[str], max_rows: int | None = None) -> list[str]:
    selected = rows if max_rows is None else rows[:max_rows]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in selected:
        lines.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return lines


def write_phase9_report(
    audit_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    key_rows: list[dict[str, str]],
    report_output: Path,
    figure_paths: dict[str, Path],
) -> None:
    tiers = Counter(row["phase8b_final_evidence_tier"] for row in audit_rows)
    classes = Counter(row["phase8b_reporting_class"] for row in audit_rows)
    high_priority = [
        row
        for row in audit_rows
        if row["gene_symbol"] in {"FTH1B", "H1F0", "RAD51", "TP53"}
    ]
    mechanism_rows = sorted(summary_rows, key=lambda row: row["mechanism"])
    lines = [
        "# Phase 9 Integrated Evidence Report",
        "",
        "Phase 9 converts the final Phase 8b evidence audit into a report and figures. It does not change evidence tiers or create new biological claims.",
        "",
        "## Primary Inputs",
        "",
        "- `results/evidence/phase8b_tier_audit.tsv`",
        "- `results/evidence/phase8b_mechanism_summary.tsv`",
        "- `results/evidence/phase8b_final_integrated_evidence.tsv`",
        "- `docs/claims_register.md`",
        "",
        "## Method Rationale",
        "",
        "The report uses deterministic categorical summaries rather than statistical testing or weighted scoring. Phase 8b rows are evidence-audit classifications, not replicate-level measurements with a sampling model. TSV-backed SVG figures are used because they are inspectable, portable, and sufficient for evidence-tier and artifact-risk summaries. Figures use a colorblind-safe palette, explicit count axes, embedded title/description metadata, and source/guardrail captions so interpretation remains traceable.",
        "",
        "## Executive Summary",
        "",
        f"- Candidates audited: {len(audit_rows)}",
        f"- Final Tier 1 rows: {tiers.get('Tier 1', 0)}",
        f"- Final Tier 2 rows: {tiers.get('Tier 2', 0)}",
        f"- Artifact/uncertain rows: {tiers.get('Artifact/uncertain', 0)}",
        f"- Plausible leads requiring validation: {classes.get('PLAUSIBLE_LEAD_REQUIRES_VALIDATION', 0)}",
        f"- Plausible leads with artifact-prone status: {classes.get('PLAUSIBLE_LEAD_ARTIFACT_PRONE', 0)}",
        "",
        "## Key Findings",
        "",
    ]
    for row in key_rows:
        lines.append(f"- **{row['question']}** {row['answer']} Guardrail: {row['interpretation_guardrail']}")
    lines.extend(
        [
            "",
            "## Figures",
            "",
            "### Evidence Tier Summary",
            "",
            f"![Evidence tier summary]({relative_for_markdown(figure_paths['tier'], report_output)})",
            "",
            "### Reporting Class Summary",
            "",
            f"![Reporting class summary]({relative_for_markdown(figure_paths['class'], report_output)})",
            "",
            "### Mechanism Evidence Matrix",
            "",
            f"![Mechanism evidence matrix]({relative_for_markdown(figure_paths['mechanism'], report_output)})",
            "",
            "### Artifact Context Summary",
            "",
            f"![Artifact context summary]({relative_for_markdown(figure_paths['artifact'], report_output)})",
            "",
            "## Mechanism Summary",
            "",
        ]
    )
    lines.extend(
        markdown_table(
            mechanism_rows,
            [
                "mechanism",
                "candidate_count",
                "tier2_count",
                "artifact_uncertain_count",
                "artifact_prone_count",
                "expression_supported_count",
                "repeat_context_caveat_count",
                "locus_review_required_count",
            ],
        )
    )
    lines.extend(
        [
            "",
            "## High-Priority Candidate Guardrails",
            "",
        ]
    )
    lines.extend(
        markdown_table(
            high_priority,
            [
                "gene_symbol",
                "mechanism",
                "phase8b_final_evidence_tier",
                "phase8b_reporting_class",
                "artifact_risk_level",
                "major_caveats",
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Interpretation Boundaries",
            "",
            "- This report can classify current rows as robust-ready, plausible, exploratory, or artifact-prone under the workflow criteria; in the current Phase 8b output, no row is robust-ready.",
            "- This report cannot infer pathway state, validated duplication, functional advantage, telomere length, organism-wide aging mechanism, causation, or human translational relevance.",
            "- Repeat context is artifact/context evidence only.",
            "- Retina expression support is tissue-specific and does not substitute for locus validation.",
            "- `TP53` remains artifact/uncertain in the current audit and should be treated as a p53-family validation problem, not a mechanism result.",
            "",
            "## Required Follow-Up",
            "",
            "- Resolve high-priority candidate loci with cross-resource support where feasible.",
            "- Add independent domain/locus validation for candidates that would otherwise be used in duplication wording.",
            "- Treat expression support as retina-specific until additional tissues and metadata are available.",
            "- Keep every biological statement traceable to the claims register and Phase 8b row-level audit.",
        ]
    )
    ensure_parent(report_output)
    report_output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_phase9_outputs(
    audit_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    report_output: Path,
    figure_manifest_output: Path,
    key_findings_output: Path,
    tier_data_output: Path,
    class_data_output: Path,
    mechanism_matrix_data_output: Path,
    artifact_context_data_output: Path,
    tier_figure_output: Path,
    class_figure_output: Path,
    mechanism_matrix_figure_output: Path,
    artifact_context_figure_output: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    if not audit_rows:
        raise ValueError("Phase 9 requires at least one Phase 8b audit row")
    if not summary_rows:
        raise ValueError("Phase 9 requires at least one Phase 8b mechanism summary row")

    tier_summary = build_tier_summary(audit_rows)
    class_summary = build_class_summary(audit_rows)
    mechanism_matrix = build_mechanism_matrix(summary_rows)
    artifact_context = build_artifact_context_rows(audit_rows, summary_rows)
    key_findings = build_key_findings(audit_rows, summary_rows)

    write_tsv(tier_data_output, tier_summary, TIER_SUMMARY_COLUMNS)
    write_tsv(class_data_output, class_summary, CLASS_SUMMARY_COLUMNS)
    write_tsv(mechanism_matrix_data_output, mechanism_matrix, MECHANISM_MATRIX_COLUMNS)
    write_tsv(artifact_context_data_output, artifact_context, ARTIFACT_CONTEXT_COLUMNS)
    write_tsv(key_findings_output, key_findings, PHASE9_KEY_FINDINGS_COLUMNS)

    write_horizontal_bar_svg(
        tier_summary,
        "evidence_tier",
        "count",
        tier_figure_output,
        "Phase 9 evidence tier summary",
        "Counts from Phase 8b final evidence tiers.",
        TIER_COLORS,
        source_label="results/evidence/phase8b_tier_audit.tsv",
        guardrail="Evidence tiers are reporting-readiness classes, not direct mechanism claims.",
    )
    write_horizontal_bar_svg(
        class_summary,
        "reporting_class",
        "count",
        class_figure_output,
        "Phase 9 reporting class summary",
        "Counts from Phase 8b reporting classes.",
        CLASS_COLORS,
        label_map=CLASS_LABELS,
        source_label="results/evidence/phase8b_tier_audit.tsv",
        guardrail="Reporting classes prioritize conservative language and do not upgrade evidence.",
    )
    write_mechanism_matrix_svg(mechanism_matrix, mechanism_matrix_figure_output)
    write_horizontal_bar_svg(
        artifact_context,
        "category",
        "count",
        artifact_context_figure_output,
        "Phase 9 artifact and context summary",
        "Counts that should harden interpretation rather than strengthen biological claims.",
        ARTIFACT_COLORS,
        label_map=ARTIFACT_LABELS,
        source_label="results/evidence/phase8b_tier_audit.tsv; results/evidence/phase8b_mechanism_summary.tsv",
        guardrail="Artifact/context counts are caution flags, not positive support for longevity mechanisms.",
    )

    figure_manifest = [
        {
            "figure_id": "FIG9-001",
            "title": "Evidence tier summary",
            "figure_type": "horizontal_bar_svg",
            "source_table": "results/evidence/phase8b_tier_audit.tsv",
            "figure_data": path_text(tier_data_output),
            "output_svg": path_text(tier_figure_output),
            "status": "GENERATED",
            "interpretation_guardrail": "Colorblind-safe horizontal bar chart with count axis; tier counts summarize reporting readiness only.",
        },
        {
            "figure_id": "FIG9-002",
            "title": "Reporting class summary",
            "figure_type": "horizontal_bar_svg",
            "source_table": "results/evidence/phase8b_tier_audit.tsv",
            "figure_data": path_text(class_data_output),
            "output_svg": path_text(class_figure_output),
            "status": "GENERATED",
            "interpretation_guardrail": "Colorblind-safe horizontal bar chart with readable class labels; reporting classes do not create biological mechanism claims.",
        },
        {
            "figure_id": "FIG9-003",
            "title": "Mechanism evidence matrix",
            "figure_type": "stacked_horizontal_bar_svg",
            "source_table": "results/evidence/phase8b_mechanism_summary.tsv",
            "figure_data": path_text(mechanism_matrix_data_output),
            "output_svg": path_text(mechanism_matrix_figure_output),
            "status": "GENERATED",
            "interpretation_guardrail": "Stacked count matrix separates evidence tiers from expression/repeat/locus context; mechanism rows require row-level caveats.",
        },
        {
            "figure_id": "FIG9-004",
            "title": "Artifact and context summary",
            "figure_type": "horizontal_bar_svg",
            "source_table": "results/evidence/phase8b_tier_audit.tsv;results/evidence/phase8b_mechanism_summary.tsv",
            "figure_data": path_text(artifact_context_data_output),
            "output_svg": path_text(artifact_context_figure_output),
            "status": "GENERATED",
            "interpretation_guardrail": "Colorblind-safe horizontal bar chart; artifact/context counts are caution flags, not positive support.",
        },
    ]
    write_tsv(figure_manifest_output, figure_manifest, PHASE9_FIGURE_MANIFEST_COLUMNS)
    write_phase9_report(
        audit_rows,
        summary_rows,
        key_findings,
        report_output,
        {
            "tier": tier_figure_output,
            "class": class_figure_output,
            "mechanism": mechanism_matrix_figure_output,
            "artifact": artifact_context_figure_output,
        },
    )
    return figure_manifest, key_findings


def run_phase9(
    config_path: Path,
    report_output: Path,
    figure_manifest_output: Path,
    key_findings_output: Path,
    tier_data_output: Path,
    class_data_output: Path,
    mechanism_matrix_data_output: Path,
    artifact_context_data_output: Path,
    tier_figure_output: Path,
    class_figure_output: Path,
    mechanism_matrix_figure_output: Path,
    artifact_context_figure_output: Path,
) -> None:
    config = read_yaml(config_path)
    if "phase9_report_generation" not in config:
        raise ValueError("config.yaml is missing phase9_report_generation")
    phase = config["phase9_report_generation"]
    audit_rows = read_tsv(Path(phase["phase8b_tier_audit"]), PHASE8B_TIER_AUDIT_COLUMNS)
    summary_rows = read_tsv(Path(phase["phase8b_mechanism_summary"]), PHASE8B_MECHANISM_SUMMARY_COLUMNS)
    generate_phase9_outputs(
        audit_rows,
        summary_rows,
        report_output,
        figure_manifest_output,
        key_findings_output,
        tier_data_output,
        class_data_output,
        mechanism_matrix_data_output,
        artifact_context_data_output,
        tier_figure_output,
        class_figure_output,
        mechanism_matrix_figure_output,
        artifact_context_figure_output,
    )
    LOGGER.info("Wrote Phase 9 report and figures from %d audited candidates", len(audit_rows))


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Phase 9 report and SVG figures from Phase 8b audit tables.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    parser.add_argument("--figure-manifest-output", type=Path, required=True)
    parser.add_argument("--key-findings-output", type=Path, required=True)
    parser.add_argument("--tier-data-output", type=Path, required=True)
    parser.add_argument("--class-data-output", type=Path, required=True)
    parser.add_argument("--mechanism-matrix-data-output", type=Path, required=True)
    parser.add_argument("--artifact-context-data-output", type=Path, required=True)
    parser.add_argument("--tier-figure-output", type=Path, required=True)
    parser.add_argument("--class-figure-output", type=Path, required=True)
    parser.add_argument("--mechanism-matrix-figure-output", type=Path, required=True)
    parser.add_argument("--artifact-context-figure-output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_phase9(
        args.config,
        args.report_output,
        args.figure_manifest_output,
        args.key_findings_output,
        args.tier_data_output,
        args.class_data_output,
        args.mechanism_matrix_data_output,
        args.artifact_context_data_output,
        args.tier_figure_output,
        args.class_figure_output,
        args.mechanism_matrix_figure_output,
        args.artifact_context_figure_output,
    )


if __name__ == "__main__":
    main()
