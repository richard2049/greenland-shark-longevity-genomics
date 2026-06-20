# Maintenance Runbook

This runbook records operational commands for refreshing validated outputs on the current native Windows laptop. It is intentionally separate from `README.md` so the project overview stays readable.

## Native Windows Policy

Use direct Python modules for periodic refreshes when possible. Broad native-Windows Snakemake dry-runs have repeatedly hung in this repository, so they are not the routine validation route here.

Before running workflow-like commands, activate the environment and work from the repository root:

```powershell
conda activate green_shark
cd D:\Documents\Greenland_Shark
```

## Focused Test Set

```powershell
python -m pytest -q `
  tests/test_phase9_package_audit.py `
  tests/test_phase9_report_generation.py `
  tests/test_phase8b_final_evidence_scoring.py `
  tests/test_phase8a_expression_evidence_integration.py `
  tests/test_phase7d_rnaseq_quantification.py `
  tests/test_phase7e_expression_hardening.py `
  tests/test_phase7c_expression_reference.py `
  tests/test_phase7b_quantification_strategy.py `
  tests/test_phase7_rnaseq.py `
  tests/test_phase6_telomere.py `
  tests/test_phase4e_locus_hardening.py `
  tests/test_phase4d_consolidation.py `
  tests/test_phase4c_locus_review.py `
  tests/test_genome_validation.py `
  --basetemp .tmp/pytest_phase_refresh
```

## Direct Refresh Through Phase 9

The commands below update current deterministic postprocessing outputs from existing validated intermediate tables. They do not run OrthoFinder, RepeatModeler, RepeatMasker, or raw-read download. Phase 7d runs Salmon only when local paired FASTQs, a Salmon executable, and a Salmon index are available; otherwise it writes explicit `NOT_RUN` status rows. Phase 7e reviews the resulting expression matrix before any expression support is allowed into Phase 8.

On native Windows, `config/config.yaml` points Phase 7d to `tools/salmon-docker.cmd`, which runs `combinelab/salmon:latest` through Docker Desktop. Existing non-empty `quant.sf` files are treated as restart checkpoints so interrupted quantification can resume without repeating successful samples.

Phase 9 is the final command in this direct refresh sequence. Phase 8b audits the Phase 8a evidence table against orthology, locus, repeat-context, telomere-readiness, and expression-support caveats. Phase 9 then renders the final conservative report and figures without changing tiers.

```powershell
python -m greenland_shark_longevity.candidate_panels `
  --input config/candidate_panels.yaml `
  --output results/validation/candidate_panel_validation.tsv

python -m greenland_shark_longevity.candidate_orthofinder `
  --candidate-panels config/candidate_panels.yaml `
  --annotation-gff "data/raw/references/SMIC_TOKYO_GENOME_2025/figshare_annotation/greenland shark annotation/complete.genomic.gff" `
  --orthogroup-gene-counts-long results/orthology/orthogroup_gene_counts_long.tsv `
  --gene-coordinates results/orthology/reference_gene_coordinates.tsv `
  --orthofinder-input-manifest results/orthology/orthofinder_input_manifest.tsv `
  --copy-number-output results/orthology/candidate_copy_number.tsv `
  --duplication-audit-output results/orthology/candidate_duplication_audit.tsv

python -m greenland_shark_longevity.evidence `
  --duplication-audit results/orthology/candidate_duplication_audit.tsv `
  --output results/evidence/phase3_integrated_evidence.tsv

python -m greenland_shark_longevity.phase4d_consolidation `
  --base-evidence results/evidence/phase3_integrated_evidence.tsv `
  --phase4c-gene-summary results/rescue/phase4c_gene_review_summary.tsv `
  --phase4c-locus-review results/rescue/phase4c_locus_review.tsv `
  --tp53-summary results/rescue/tp53_targeted_forward_search_summary.tsv `
  --interpretation-output results/evidence/phase4d_candidate_interpretation.tsv `
  --integrated-output results/evidence/phase4d_integrated_evidence.tsv

python -m greenland_shark_longevity.phase4e_locus_hardening `
  --phase4c-locus-review results/rescue/phase4c_locus_review.tsv `
  --domain-integrity results/rescue/phase4c_rescue_domain_integrity.tsv `
  --tp53-alignment-hits results/rescue/tp53_targeted_forward_alignment_hits.tsv `
  --tp53-target-regions results/rescue/tp53_forward_hit_target_regions.tsv `
  --annotation-gff "data/raw/references/SMIC_TOKYO_GENOME_2025/figshare_annotation/greenland shark annotation/complete.genomic.gff" `
  --locus-output results/rescue/phase4e_locus_manual_review.tsv `
  --summary-output results/rescue/phase4e_gene_hardened_summary.tsv `
  --evidence-output results/evidence/phase4e_hardened_evidence.tsv

python -m greenland_shark_longevity.phase5_repeat_context `
  --config config/config.yaml `
  --candidate-loci results/rescue/phase4e_locus_manual_review.tsv `
  --figshare-inventory data/metadata/figshare_annotation_inventory.tsv `
  --resource-status-output data/metadata/phase5_repeat_resource_status.tsv `
  --repeat-features-output results/repeats/phase5_repeat_features.tsv `
  --locus-context-output results/repeats/phase5_candidate_locus_repeat_context.tsv `
  --gene-summary-output results/repeats/phase5_gene_repeat_context_summary.tsv

python -m greenland_shark_longevity.phase5c_repeat_qc `
  --config config/config.yaml `
  --candidate-loci results/rescue/phase4e_locus_manual_review.tsv `
  --phase5-locus-context results/repeats/phase5_candidate_locus_repeat_context.tsv `
  --phase5b-inventory data/metadata/phase5b_repeat_output_inventory.tsv `
  --integrity-output results/repeats/phase5c_repeatmasker_integrity.tsv `
  --locus-qc-output results/repeats/phase5c_locus_repeat_qc.tsv `
  --gene-qc-output results/repeats/phase5c_gene_repeat_qc_summary.tsv `
  --report-output reports/generated/phase5c_repeat_context_qc_report.md

python -m greenland_shark_longevity.phase5_evidence_integration `
  --base-evidence results/evidence/phase4d_integrated_evidence.tsv `
  --phase4e-evidence results/evidence/phase4e_hardened_evidence.tsv `
  --phase5-gene-summary results/repeats/phase5_gene_repeat_context_summary.tsv `
  --phase5c-gene-qc results/repeats/phase5c_gene_repeat_qc_summary.tsv `
  --phase5-evidence-output results/evidence/phase5_repeat_context_evidence.tsv `
  --integrated-output results/evidence/integrated_evidence.tsv `
  --report-output reports/generated/phase5_repeat_context_report.md

python -m greenland_shark_longevity.phase6_telomere `
  --config config/config.yaml `
  --motif-scan-output results/telomere/phase6_telomeric_motif_scan.tsv `
  --enrichment-output results/telomere/phase6_scaffold_end_enrichment.tsv `
  --gene-audit-output results/telomere/phase6_telomere_gene_audit.tsv `
  --report-output reports/generated/phase6_telomere_report.md

python -m greenland_shark_longevity.phase7_rnaseq `
  --config config/config.yaml `
  --manifest-output data/metadata/rnaseq_manifest.tsv `
  --readiness-output results/rnaseq/phase7_rnaseq_readiness.tsv `
  --expression-plan-output results/rnaseq/phase7_candidate_expression_plan.tsv `
  --report-output reports/generated/phase7_rnaseq_readiness_report.md

python -m greenland_shark_longevity.phase7b_quantification_strategy `
  --config config/config.yaml `
  --reference-strategy-output results/rnaseq/phase7b_reference_quantification_strategy.tsv `
  --candidate-quant-map-output results/rnaseq/phase7b_candidate_quantification_map.tsv `
  --run-plan-output results/rnaseq/phase7b_quantification_run_plan.tsv `
  --report-output reports/generated/phase7b_quantification_strategy_report.md

python -m greenland_shark_longevity.phase7c_expression_reference `
  --config config/config.yaml `
  --transcript-fasta-output data/interim/rnaseq/phase7c/smic_tokyo.annotation_transcripts.fna `
  --tx2gene-output data/interim/rnaseq/phase7c/smic_tokyo.tx2gene.tsv `
  --reference-qc-output results/rnaseq/phase7c_expression_reference_qc.tsv `
  --candidate-validation-output results/rnaseq/phase7c_candidate_reference_validation.tsv `
  --report-output reports/generated/phase7c_expression_reference_report.md

python -m greenland_shark_longevity.phase7d_rnaseq_quantification `
  --config config/config.yaml `
  --raw-read-intake-output results/rnaseq/phase7d_raw_read_intake.tsv `
  --fastq-qc-output results/rnaseq/phase7d_fastq_qc.tsv `
  --salmon-preflight-output results/rnaseq/phase7d_salmon_preflight.tsv `
  --salmon-quant-summary-output results/rnaseq/phase7d_salmon_quant_summary.tsv `
  --candidate-expression-matrix-output results/rnaseq/phase7d_candidate_expression_matrix.tsv `
  --report-output reports/generated/phase7d_rnaseq_quantification_report.md

python -m greenland_shark_longevity.phase7e_expression_hardening `
  --config config/config.yaml `
  --run-qc-output results/rnaseq/phase7e_run_qc_review.tsv `
  --candidate-output results/rnaseq/phase7e_candidate_expression_hardened.tsv `
  --parameter-review-output results/rnaseq/phase7e_parameter_review.tsv `
  --report-output reports/generated/phase7e_expression_hardening_report.md

python -m greenland_shark_longevity.phase8a_expression_evidence_integration `
  --config config/config.yaml `
  --expression-evidence-output results/evidence/phase8a_expression_support_evidence.tsv `
  --integration-audit-output results/evidence/phase8a_expression_integration_audit.tsv `
  --integrated-output results/evidence/phase8a_integrated_evidence.tsv `
  --report-output reports/generated/phase8a_expression_evidence_report.md

python -m greenland_shark_longevity.phase8b_final_evidence_scoring `
  --config config/config.yaml `
  --final-evidence-output results/evidence/phase8b_final_integrated_evidence.tsv `
  --tier-audit-output results/evidence/phase8b_tier_audit.tsv `
  --mechanism-summary-output results/evidence/phase8b_mechanism_summary.tsv `
  --report-output reports/generated/phase8b_final_evidence_report.md

python -m greenland_shark_longevity.phase9_report_generation `
  --config config/config.yaml `
  --report-output reports/final/greenland_shark_longevity_phase9_report.md `
  --figure-manifest-output results/reporting/phase9_figure_manifest.tsv `
  --key-findings-output results/reporting/phase9_key_findings.tsv `
  --tier-data-output reports/figures/data/phase9_evidence_tier_summary.tsv `
  --class-data-output reports/figures/data/phase9_reporting_class_summary.tsv `
  --mechanism-matrix-data-output reports/figures/data/phase9_mechanism_evidence_matrix.tsv `
  --artifact-context-data-output reports/figures/data/phase9_artifact_context_summary.tsv `
  --tier-figure-output reports/figures/phase9_evidence_tier_summary.svg `
  --class-figure-output reports/figures/phase9_reporting_class_summary.svg `
  --mechanism-matrix-figure-output reports/figures/phase9_mechanism_evidence_matrix.svg `
  --artifact-context-figure-output reports/figures/phase9_artifact_context_summary.svg

python -m greenland_shark_longevity.phase9_package_audit `
  --config config/config.yaml `
  --package-audit-output results/reporting/phase9_report_package_audit.tsv `
  --release-readiness-output results/reporting/phase9_public_repository_readiness.tsv `
  --report-output reports/final/phase9_report_package_audit.md
```

## Managed Snakemake Target

Use this only as a controlled single command, preferably in WSL/Linux, Docker-backed environments, or when direct modules are insufficient:

```powershell
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=phase9_report_package_audit
```

Earlier stopping points remain available:

- `workflow_mode=phase5_repeat_context`
- `workflow_mode=phase5c_repeat_qc`
- `workflow_mode=phase6_telomere`
- `workflow_mode=phase7_rnaseq_readiness`
- `workflow_mode=phase7b_quantification_strategy`
- `workflow_mode=phase7c_expression_reference`
- `workflow_mode=phase7d_rnaseq_quantification`
- `workflow_mode=phase7e_expression_hardening`
- `workflow_mode=phase8a_expression_integration`
- `workflow_mode=phase8b_final_evidence_scoring`
- `workflow_mode=phase9_report_generation`

## Interruption Recovery

If a native-Windows workflow command is interrupted or appears frozen, inspect exact process command lines before stopping anything:

```powershell
Get-CimInstance Win32_Process |
  Where-Object { $_.Name -match 'python|snakemake' -or $_.CommandLine -match 'snakemake|greenland_shark_longevity' } |
  Select-Object ProcessId,Name,CommandLine
```

Stop only verified orphaned Snakemake/Python PIDs that match the interrupted repository command.

## WSL/Linux OrthoFinder Route

If OrthoFinder is available inside WSL/Linux, the integrated Snakemake mode can run OrthoFinder and downstream parsing:

```bash
cd /mnt/d/Documents/Greenland_Shark
conda env create -f environment.orthology.yml
conda activate green_shark_orthology
orthofinder -h
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=orthology
```

On native Windows, the conda environment does not provide OrthoFinder because the Bioconda package depends on external search tools such as `blast`, which are not available for this `win-64` solve. The practical Windows route is Docker OrthoFinder followed by repository parsing and postprocessing. See `docs/orthology_setup.md`.
