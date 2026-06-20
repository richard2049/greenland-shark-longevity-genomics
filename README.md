# Greenland Shark Longevity Genomics

Reproducible workflow for evaluating public genomic and transcriptomic evidence related to extreme longevity in the Greenland shark (*Somniosus microcephalus*).

This repository is designed for careful interpretation. It does not treat public annotations as sufficient evidence for a longevity mechanism. Instead, it asks which candidate signals are supported by current resources, which are plausible but incomplete, and which remain too artifact-prone to claim biologically.

## Why This Exists

Greenland sharks are among the longest-lived vertebrates, but public genomic resources are still uneven: assemblies differ, annotations are incomplete, gene families can be hard to resolve, and transcriptomic data are limited in tissue scope.

This project turns those limitations into explicit audit fields rather than hiding them. The workflow tracks provenance, tests candidate panels, compares orthogroups, checks domains and loci, records repeat context, reviews retina RNA-seq support, and produces reports that keep interpretation guardrails visible.

## Current Scope

The repository currently includes a complete v0.1 workflow through Phase 9:

| Area | Current implementation |
| --- | --- |
| Metadata and provenance | Registered public Greenland shark genome/RNA-seq resources and comparator proteomes |
| Resource QC | Assembly, protein, annotation, file-inventory, and metadata checks |
| Candidate panels | DNA repair, p53, chromatin, telomere, ferroptosis/iron, antioxidant response, mitochondria, immune regulation, proteostasis, retina |
| Orthology and copy number | OrthoFinder postprocessing, candidate copy-number tables, duplication-audit tables |
| Annotation rescue | Targeted high-priority candidate review for `H1F0`, `FTH1B`, `RAD51`, and `TP53` |
| Domains and loci | HMMER/Pfam domain checks, miniprot-assisted locus review, manual hardening tables |
| Repeat context | Candidate-locus repeat overlap and repeat-context QC, reported only as artifact/context evidence |
| Telomere layer | Canonical motif scan and telomere-gene readiness checks, not telomere-length inference |
| RNA-seq support | Retina-only metadata, reference construction, Salmon quantification support, expression hardening |
| Evidence scoring | Final Phase 8b audit and Phase 9 report/figures |

## Current Evidence Snapshot

The current Phase 9 report audits 41 candidate rows:

- Tier 1: 0 rows
- Tier 2: 31 rows
- Artifact/uncertain: 10 rows
- Artifact-prone plausible leads: `FTH1B`, `H1F0`, `RAD51`
- `TP53`: remains `Artifact/uncertain` because the current p53-family signal is not resolved enough for a gene-state or mechanism claim
- Retina expression support: tissue-specific context only, not differential expression or pathway activity

These are workflow audit results, not standalone biological conclusions. See `reports/final/greenland_shark_longevity_phase9_report.md` and `results/evidence/phase8b_tier_audit.tsv`.

## What This Repository Does Not Claim

- No gene loss from annotation absence.
- No pathway activation or inactivation without an appropriate design.
- No validated duplication from OrthoFinder counts alone.
- No telomere length or telomerase activity from motif scans.
- No repeat-mediated longevity mechanism from repeat overlap.
- No whole-organism aging interpretation from retina-only RNA-seq.
- No human translational relevance outside explicitly speculative future work.

## Quickstart

Create the environment:

```bash
conda env create -f environment.yml
conda activate green_shark
```

Run the lightweight demo workflow:

```bash
snakemake --snakefile workflow/Snakefile --cores 1
```

Run metadata/provenance only:

```bash
snakemake --snakefile workflow/Snakefile --cores 1 --config workflow_mode=metadata_only
```

Run the final report-package audit:

```powershell
python -m greenland_shark_longevity.phase9_package_audit `
  --config config/config.yaml `
  --package-audit-output results/reporting/phase9_report_package_audit.tsv `
  --release-readiness-output results/reporting/phase9_public_repository_readiness.tsv `
  --report-output reports/final/phase9_report_package_audit.md
```

Native-Windows note: broad Snakemake dry-runs have been unreliable on this laptop. For maintenance, prefer direct Python module commands from `docs/runbook.md`; use Docker, WSL/Linux, or external compute for Linux-first tools such as OrthoFinder, miniprot, Salmon, RepeatModeler, and RepeatMasker.

## Key Outputs

Small final release outputs are tracked in Git:

- `results/evidence/phase8b_final_integrated_evidence.tsv`
- `results/evidence/phase8b_tier_audit.tsv`
- `results/evidence/phase8b_mechanism_summary.tsv`
- `results/reporting/phase9_key_findings.tsv`
- `results/reporting/phase9_figure_manifest.tsv`
- `reports/final/greenland_shark_longevity_phase9_report.md`
- `reports/final/phase9_report_package_audit.md`
- `reports/figures/*.svg`
- `reports/figures/data/*.tsv`

Large raw data, intermediates, logs, Salmon outputs, repeat outputs, OrthoFinder working folders, and workflow caches are intentionally ignored.

## Repository Layout

```text
config/      Workflow configuration, resources, and candidate panels
data/demo/   Tiny non-biological fixtures for runnable tests and demo mode
data/metadata/
             Public-resource manifests and file inventories
docs/        Study design, claims register, runbook, setup notes, release checklist
reports/     Final report, figure data, and rendered SVG summaries
results/     Small final evidence/reporting TSV release bundle
src/         Python package and command-line modules
tests/       Pytest coverage for schemas, parsers, scoring, and report generation
tools/       Small local wrappers for Docker-based tools
workflow/    Snakemake entry point
```

## Documentation Guide

- `docs/study_design.md`: scientific questions, phase model, evidence tiers, limitations
- `docs/claims_register.md`: claim-level traceability and caveats
- `docs/data_manifest.md`: registered public resources and usage notes
- `docs/runbook.md`: direct refresh commands for deterministic postprocessing
- `docs/orthology_setup.md`: OrthoFinder setup and interpretation
- `docs/domain_validation.md`: HMMER/Pfam domain-validation notes
- `docs/annotation_rescue.md`: targeted rescue and locus-validation notes
- `docs/repeat_annotation_setup.md`: optional RepeatModeler/RepeatMasker route
- `docs/release_checklist.md`: public-release checklist
- `docs/maintainability_review.md`: technical review and known debt

## Data Policy

This repository is designed to be public without committing large biological data.

Tracked:

- Code
- Config
- Documentation
- Demo fixtures
- Metadata manifests
- Small final evidence/reporting TSVs
- Final report and figure files

Ignored:

- Raw genome, protein, RNA-seq, and repeat files
- FASTQ/BAM/CRAM/SAM and other large sequence files
- `data/raw/`
- `data/interim/`
- `logs/`
- `.snakemake/`
- Heavy tool outputs and local indexes

## Development Checks

Run the full test suite:

```bash
python -m pytest -q
```

Run focused Phase 9 checks:

```bash
python -m pytest -q tests/test_phase9_report_generation.py tests/test_phase9_package_audit.py
```

Run syntax/import sanity:

```bash
python -m compileall -q src/greenland_shark_longevity
```

## Citation And License

The repository code and documentation are released under the MIT License; see `LICENSE`.

If you use or adapt this workflow, cite it using `CITATION.cff`. Public biological datasets referenced by the workflow retain their original database, publication, and provider terms.
