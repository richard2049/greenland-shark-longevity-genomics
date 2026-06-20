# Maintainability Review

This review records the June 2026 cleanup pass over the repository structure, documentation, workflow entry points, and Phase 7 implementation. It is a software-quality note, not a biological result.

## Scope Reviewed

- `README.md`
- `docs/*.md`
- `config/config.yaml`
- `workflow/Snakefile`
- `src/greenland_shark_longevity/*.py`
- `tests/test_phase7*.py`
- Generated cache directories visible in the workspace

Large raw data, result tables, provenance files, and scientific intermediate outputs were not deleted. They are part of the reproducibility record unless explicitly superseded by a documented workflow step.

## Changes Made

- Moved the long periodic refresh command block out of `README.md` into `docs/runbook.md`.
- Reordered the Phase 5b README section so it appears with the Phase 5 repeat-analysis material rather than after Phase 7.
- Kept README focused on project purpose, core commands, phase summaries, expected outputs, and pointers to detailed docs.
- Removed generated Python and pytest cache directories from the workspace after verifying their paths were inside the repository.
- Kept Phase 7c as a standalone deterministic module rather than embedding transcript-reference construction into notebooks or ad hoc shell commands.
- Centralized repeated Phase 7 helper logic in `src/greenland_shark_longevity/utils.py`: `NOT_ASSESSED`, gzip-aware text opening, text cleanup, delimited-value splitting, and TSV-safe value joining.
- Updated Phase 7a, Phase 7b, and Phase 7c modules to use the shared helpers while preserving output schemas.
- Added Phase 7d as a guarded, retina-only raw-read intake/QC and Salmon quantification module with explicit `NOT_RUN` statuses when local reads, Salmon, or an index are unavailable.

## Script Review Findings

- The largest modules are Phase 4/5 validation and repeat-context scripts. They are long, but they are phase-specific, tested, and tied to stable output schemas. They should not be merged or rewritten without a dedicated refactor plan and golden-output tests.
- Phase 7 scripts had duplicated small helper functions. Those were low-risk to consolidate because they do not encode scientific decisions.
- `workflow/Snakefile` uses repeated mode-condition sets. This is readable enough for the current MVP, but it is a future refactor candidate if additional phases are added.
- The source tree had generated Python and pytest cache directories. These are ignored by `.gitignore` and were removed from the workspace.
- Source strings that include restricted biological terms such as `activated`, `inactivated`, `absent`, `causes`, and `proves` are guardrail strings or tests that prevent overclaiming.

## Review Decisions

- No scientific result files were deleted. Even imperfect outputs such as repeat-context QC tables and warning-bearing logs remain useful provenance.
- No broad refactor was applied to older Phase 4/5 modules. They are larger than ideal, but they encode validated phase behavior and have tests. Refactoring them now would add risk without improving the immediate Phase 7 objective.
- The current `README.md` is still substantial, but it now avoids the largest operational command block. Further reduction should happen only after creating stable per-phase docs for Phase 4, Phase 5, Phase 6, and Phase 7.
- Generated transcript FASTA and `tx2gene` files are kept under `data/interim/` because they are reproducible intermediates needed for later RNA-seq quantification.

## Remaining Technical Debt

- `workflow/Snakefile` has long mode-condition sets. A future cleanup could centralize phase dependencies in a small Python helper or Snakemake module, but this is not urgent while native-Windows Snakemake remains secondary.
- `config/config.yaml` is long because it stores provenance and phase settings together. A future split into `config/resources.yaml`, `config/phases.yaml`, and `config/runtime.yaml` may improve readability, but should be done only with schema tests.
- Phase 4 and Phase 5 modules are large. Future refactoring should preserve exact output schemas and be driven by tests, not style alone.
- Future Phase 7 work should add full read-QC reports and locus-aware genome-aligned validation only after local FASTQ provenance, Salmon index provenance, and candidate quantification outputs are stable.

## Scientific Guardrail

This cleanup does not upgrade evidence tiers or add biological claims. Documentation and code wording should continue to distinguish:

- resource availability,
- input readiness,
- candidate reference mapping,
- expression quantification,
- expression interpretation,
- and longevity-mechanism evidence.

## Final Report And Release-Readiness Review

This follow-up review checked the Phase 9 report package after the final evidence/reporting TSV release-bundle policy was added.

Reviewed files and surfaces:

- `README.md`
- `.gitignore`
- `docs/release_checklist.md`
- `docs/claims_register.md`
- `reports/final/greenland_shark_longevity_phase9_report.md`
- `reports/final/phase9_report_package_audit.md`
- `results/evidence/phase8b_final_integrated_evidence.tsv`
- `results/evidence/phase8b_tier_audit.tsv`
- `results/evidence/phase8b_mechanism_summary.tsv`
- `results/reporting/phase9_figure_manifest.tsv`
- `results/reporting/phase9_key_findings.tsv`
- `src/greenland_shark_longevity/phase9_report_generation.py`
- `src/greenland_shark_longevity/phase9_package_audit.py`

Decisions from this review:

- Kept the Phase 8b and Phase 9 filenames because the phase prefix is useful provenance, the names are stable across config/tests/docs, and renaming them immediately before public packaging would weaken traceability.
- Kept raw data, intermediates, Salmon outputs, repeat outputs, logs, and broad workflow products ignored. Only the small final evidence/reporting TSV release bundle is unignored by exact `.gitignore` exceptions.
- Kept the large Phase 4/5/7 scripts as phase-specific modules. Several exceed the ideal size for long-term maintenance, but they are schema-driven and tested. Splitting them should be a dedicated refactor with golden-output tests, not a cosmetic pre-release edit.
- Clarified Phase 9 wording so the final report says no row is currently robust-ready instead of implying that robust findings are present.
- Confirmed that the final report keeps repeat context, retina expression support, and TP53/p53-family uncertainty inside explicit interpretation boundaries.

Remaining release-readiness cautions:

- The repository folder is not currently initialized as a Git repository, so staged-file checks cannot yet be performed with `git status`.
- Generated Python caches and egg-info directories can reappear during tests or editable installs. They are ignored by `.gitignore` and should not be part of a public release bundle.
- The current SVG figure stack is acceptable for repository evidence-audit reporting. A future manuscript-oriented figure package could add optional matplotlib, seaborn, ggplot2, or R export layers that read the same figure-data TSVs.
