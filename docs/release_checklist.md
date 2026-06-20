# Public Release Checklist

This checklist prepares the Greenland shark longevity genomics repository for public sharing. It is a repository-quality checklist, not a biological validation result.

## Scientific Scope

- [ ] Confirm the final report still states that no current Phase 8b row reaches Tier 1.
- [ ] Confirm `FTH1B`, `H1F0`, and `RAD51` remain plausible but artifact-prone unless new validation changes the evidence tables.
- [ ] Confirm `TP53` remains `Artifact/uncertain` unless a future validated locus/paralog analysis changes the evidence tables.
- [ ] Confirm repeat context is described only as artifact/context evidence.
- [ ] Confirm retina expression support is described only as tissue-specific support.
- [ ] Confirm no report text claims pathway state, validated duplication, telomere length, causation, functional advantage, or human translational relevance.

## Traceability

- [ ] Regenerate Phase 8b and Phase 9 outputs from the runbook if upstream tables changed.
- [ ] Run the Phase 9 package audit:

```powershell
python -m greenland_shark_longevity.phase9_package_audit `
  --config config/config.yaml `
  --package-audit-output results/reporting/phase9_report_package_audit.tsv `
  --release-readiness-output results/reporting/phase9_public_repository_readiness.tsv `
  --report-output reports/final/phase9_report_package_audit.md
```

- [ ] Confirm `results/reporting/phase9_report_package_audit.tsv` has no `FAIL` rows.
- [ ] Confirm `results/reporting/phase9_public_repository_readiness.tsv` has no unresolved `WARN` rows, or that any warning is explicitly documented.
- [ ] Confirm the final report cites `docs/claims_register.md`.

## Release Bundle

The repository intentionally keeps raw data, intermediates, Salmon outputs, repeat outputs, logs, and broad result directories ignored. Only the small final evidence/reporting TSV release bundle is unignored by exact `.gitignore` exceptions.

Release-bundle TSVs:

- [ ] `results/evidence/phase8b_final_integrated_evidence.tsv`
- [ ] `results/evidence/phase8b_tier_audit.tsv`
- [ ] `results/evidence/phase8b_mechanism_summary.tsv`
- [ ] `results/reporting/phase9_figure_manifest.tsv`
- [ ] `results/reporting/phase9_key_findings.tsv`
- [ ] `results/reporting/phase9_report_package_audit.tsv`
- [ ] `results/reporting/phase9_public_repository_readiness.tsv`

Report/figure release artifacts:

- [ ] `reports/final/greenland_shark_longevity_phase9_report.md`
- [ ] `reports/final/phase9_report_package_audit.md`
- [ ] `reports/figures/phase9_evidence_tier_summary.svg`
- [ ] `reports/figures/phase9_reporting_class_summary.svg`
- [ ] `reports/figures/phase9_mechanism_evidence_matrix.svg`
- [ ] `reports/figures/phase9_artifact_context_summary.svg`
- [ ] `reports/figures/data/phase9_evidence_tier_summary.tsv`
- [ ] `reports/figures/data/phase9_reporting_class_summary.tsv`
- [ ] `reports/figures/data/phase9_mechanism_evidence_matrix.tsv`
- [ ] `reports/figures/data/phase9_artifact_context_summary.tsv`

## Validation

- [ ] Run focused Phase 9 tests:

```powershell
python -m pytest -q `
  tests/test_phase9_package_audit.py `
  tests/test_phase9_report_generation.py `
  --basetemp .tmp/pytest_phase9_release
```

- [ ] Run the full test suite:

```powershell
python -m pytest -q --basetemp .tmp/pytest_release_full
```

- [ ] Do not use a broad native-Windows Snakemake dry-run as the default release validation route. Use direct module commands from `docs/runbook.md`, or run Snakemake in WSL/Linux/Docker if orchestration validation is required.

## Data Hygiene

- [ ] Confirm `data/raw/`, `data/interim/`, `logs/`, `.snakemake/`, and large sequencing files remain ignored.
- [ ] Confirm no FASTQ, BAM/CRAM/SAM, genome FASTA, Salmon index, repeat-model database, or large Docker/workflow artifact is staged for release.
- [ ] Confirm external resources retain accession/version/provenance records in metadata/config/docs.

## Final Review

- [ ] Read the final report once as a skeptical reviewer and check every biological statement against a supporting TSV or claim-register row.
- [ ] Confirm the README quickstart and runbook commands are synchronized.
- [ ] Confirm the repository does not require private data, hardcoded local absolute paths, or undocumented manual files for the default/demo workflow.
