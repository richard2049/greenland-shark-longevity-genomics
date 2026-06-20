# Greenland Shark Repo-Specific Agent Instructions

These instructions supplement the global computational biology and bioinformatics working agreements. They are specific to this repository and to the observed behavior of the workflow on the current native Windows laptop.

## Process-Control Rules

- Do not run more than one Snakemake command at a time in this repository.
- Do not run Snakemake commands through parallel tool execution.
- Do not run a Snakemake dry-run and a Snakemake execution command concurrently in this repository.
- Avoid broad Snakemake validation when a direct Python module command can validate the same frozen output with narrower scope and clearer logs.
- If a Snakemake command has previously frozen, do not rerun the same broad command without changing the execution strategy.
- Prefer direct, deterministic Python module commands for frozen postprocessing and evidence consolidation.
- Prefer visible execution from the already activated `green_shark` environment over buffered wrappers such as `conda run`, unless a wrapper is required for reproducibility.
- Use explicit timeouts for commands that may hang.
- Treat raw data and schema/test-validated intermediate outputs as checkpoints. Resume from checkpoints instead of rerunning heavy or unstable workflow steps when possible.

## Native Windows Snakemake Policy

This repository has shown unstable behavior with Snakemake commands on native Windows, including dry-runs that can hang, leave orphaned `snakemake.exe`/`python.exe` child processes, or trigger Windows permission errors during process inspection. On this platform:

- Treat native-Windows Snakemake dry-runs as non-default and potentially unsafe for productivity in this repository.
- Do not use Snakemake dry-runs as routine validation after code edits.
- Do not run broad Phase 4, Phase 5, Phase 6, Phase 7, or later Snakemake dry-runs on native Windows unless the user explicitly asks for that exact check or no direct validation route exists.
- Use Snakemake for controlled single-command phase execution only when it is the clearest reproducible route for the current phase and previous process state has been checked.
- Do not use broad or parallel Snakemake as a routine validation tool on native Windows for this repository.
- When validating Phase 4d and later deterministic postprocessing phases, use direct Python module commands rather than broad workflow execution.
- If full workflow orchestration is needed, prefer WSL/Linux, Docker, or an external compute environment where the relevant tools are better supported and where command logs can be captured cleanly.

## Preferred Native-Windows Validation Pattern

When implementing or refreshing a deterministic phase on this laptop, use this order instead of a Snakemake dry-run:

1. Run the smallest relevant pytest files with an explicit `--basetemp` under `.tmp/`.
2. Run the direct `python -m greenland_shark_longevity.<module>` command for the phase.
3. Check the expected output files with `Get-Item` and, when useful, `Get-Content -TotalCount`.
4. Run `python -m compileall -q src\greenland_shark_longevity` for syntax/import sanity.
5. Record in the response that broad Snakemake dry-run was intentionally skipped on native Windows because direct module validation is the safer repository policy.

Use Snakemake on native Windows only for a bounded execution target when the command itself is part of the requested workflow and the direct module route is insufficient. If a Snakemake check is unavoidable:

- Run it as a single foreground command, never through parallel tool execution.
- Use `--cores 1` for validation checks unless there is a specific reason to use more.
- Use an explicit timeout.
- Check for existing `snakemake.exe`, `python.exe`, or phase-specific tool processes before starting.
- If the command hangs or is interrupted, inspect exact command lines before stopping only the verified orphaned PIDs.
- Prefer WSL/Linux or Docker for DAG validation and full orchestration if native Windows repeats the same failure mode.

## Recovery After A Hang Or Interruption

After a user interruption, timeout, or frozen command:

1. First try to check exact Snakemake/Python command lines:

   ```powershell
   Get-CimInstance Win32_Process |
     Where-Object { $_.Name -match 'python|snakemake' -or $_.CommandLine -match 'snakemake|greenland_shark_longevity' } |
     Select-Object ProcessId,Name,CommandLine
   ```

2. If command-line inspection is blocked by Windows permissions, request permission for a read-only process check rather than guessing.

3. If command-line inspection is not available, use CPU/process listing only as a weaker fallback:

   ```powershell
   Get-Process -Name python -ErrorAction SilentlyContinue | Select-Object Id,CPU
   ```

4. Stop only verified orphaned processes related to the current task:

   ```powershell
   Stop-Process -Id <PID> -Force
   ```

5. Do not use blanket process-kill commands.
6. Do not stop a process unless its PID, command line or command type, and likely relationship to the current task are clear.
7. Check expected outputs directly before rerunning any workflow command.
8. Check `.snakemake\locks` only after verifying no active Snakemake process remains.
9. Document which PIDs were stopped and why.

## Preferred Manual Phase 4d Route

For frozen Phase 4 consolidation, use this direct route from PowerShell with `green_shark` active:

```powershell
conda activate green_shark
cd D:\Documents\Greenland_Shark

python -m pytest -q `
  tests/test_phase4d_consolidation.py `
  tests/test_phase4c_locus_review.py `
  tests/test_genome_validation.py `
  --basetemp .tmp/pytest_phase4d_related

python -m greenland_shark_longevity.evidence `
  --duplication-audit results/orthology/candidate_duplication_audit.tsv `
  --output results/evidence/phase3_integrated_evidence.tsv

python -m greenland_shark_longevity.phase4d_consolidation `
  --base-evidence results/evidence/phase3_integrated_evidence.tsv `
  --phase4c-gene-summary results/rescue/phase4c_gene_review_summary.tsv `
  --phase4c-locus-review results/rescue/phase4c_locus_review.tsv `
  --tp53-summary results/rescue/tp53_targeted_forward_search_summary.tsv `
  --interpretation-output results/evidence/phase4d_candidate_interpretation.tsv `
  --integrated-output results/evidence/integrated_evidence.tsv
```

Expected Phase 4d outputs:

- `results/evidence/phase3_integrated_evidence.tsv`
- `results/evidence/phase4d_candidate_interpretation.tsv`
- `results/evidence/integrated_evidence.tsv`

## Scientific Interpretation Guardrail

A frozen/manual consolidation step must not introduce new biological claims. It may only:

- preserve pre-consolidation evidence,
- consolidate existing validation tables,
- update evidence tiers conservatively,
- record caveats and required validation,
- and separate resource-quality observations from biological interpretation.

For the current Phase 4d state:

- `H1F0`, `FTH1B`, and `RAD51` are Tier 2 candidate-locus findings that require manual locus review and cross-resource support before duplication or functional language is justified.
- `FTH1B` remains high artifact risk because ferritin-family paralogs, isoforms, fragments, and local genomic context need manual resolution.
- `TP53` remains `Artifact/uncertain`; the targeted p53-family alignment is not a gene-state, function, or mechanism claim.

## Phase Discipline

- Keep work aligned with the repository phase model:
  - Phase 0: metadata and provenance.
  - Phase 1: assembly and annotation QC.
  - Phase 2: candidate panels.
  - Phase 3: orthology and copy-number analysis.
  - Phase 4: targeted annotation rescue and validation.
  - Phase 5: repeats and transposable elements.
  - Phase 6: telomere-related analysis.
  - Phase 7: RNA-seq support.
  - Phase 8: integrated evidence scoring.
  - Phase 9: reports and figures.
- Do not jump to heavy downstream analyses when prerequisite manifests, QC tables, schemas, and validation outputs are missing.
- Keep discovery, validation, interpretation, and reporting as separate steps with separate outputs.
- If an output is generated from demo or artificial data, label it as `DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE` or equivalent.
- The default repository state should remain runnable without large downloads or heavy tools.

## AI-Assisted Research Rules

- AI-generated code, prose, tables, or interpretations must be treated as draft material until validated against repository files or primary sources.
- Do not invent accessions, URLs, retrieval dates, checksums, filenames, sample metadata, gene functions, tool versions, database versions, or biological conclusions.
- If a fact cannot be verified locally or from a primary source, write `TODO`, `NOT_ASSESSED`, `ANNOTATION_UNCERTAINTY`, or `REQUIRES_VALIDATION`.
- Biological statements must point to machine-readable supporting files whenever possible.
- Do not use AI-generated summaries of papers as evidence unless the relevant claim has been checked against the paper, supplement, database record, or generated result table.
- Do not let fluent language upgrade weak evidence. If the data support only a candidate, uncertainty, or artifact-prone signal, say that directly.
- When changing scientific logic, update tests and documentation in the same task where practical.

## Provenance And FAIR Data Rules

- Every external resource used in analysis should have, when available:
  - accession or stable identifier,
  - version or release,
  - source URL,
  - retrieval date,
  - expected local path,
  - file type,
  - checksum or file size,
  - usage notes,
  - associated publication or resource.
- Keep raw files under ignored raw-data directories and treat them as immutable.
- Do not edit raw downloaded files in place. Create derived files under `data/interim/`, `results/`, or another documented derived-output directory.
- Do not commit large raw data, temporary indexes, workflow caches, or container outputs unless explicitly intended and documented.
- If a public resource changes or is re-downloaded, preserve the old provenance record or explicitly record the replacement.
- Prefer stable, tabular, machine-readable outputs with explicit schemas over prose-only results.

## Reproducibility And Runtime Rules

- Record tool versions, database versions, command-line parameters, thread counts, random seeds, and container image tags where relevant.
- Containerized tools are acceptable and professional when native installation is unreliable, especially for bioinformatics tools with Linux-first packaging. Record the image name/tag and command.
- Do not mix outputs from different tool versions, database versions, or parameters under the same filename unless the output is intentionally regenerated and documented.
- Prefer deterministic postprocessing scripts over manual spreadsheet edits.
- If manual inspection is required, record what was inspected, where it came from, and what was not assessed.
- Use logs for long-running external tools and retain stderr when it contains useful provenance or warnings.

## Orthology, Copy Number, And Gene-Family Rules

- Do not treat a gene symbol as sufficient evidence of orthology.
- Do not equate "not annotated" with gene absence.
- Do not equate copy-number expansion with functional advantage.
- Do not infer gene loss or gene-state changes without direct sequence-level evidence and manual review.
- Candidate duplication claims require, as feasible:
  - orthology support,
  - isoform filtering,
  - separable genomic loci,
  - protein/domain integrity,
  - coordinate support,
  - cross-resource consistency,
  - and local genomic-context or synteny review.
- Raw OrthoFinder counts are evidence inputs, not final copy-number claims.
- For complex families such as histones, ferritins, RAD51 paralogs, immune genes, and p53-family genes, assume paralogy and annotation ambiguity until specifically resolved.
- If repeat annotations are available, inspect repeat-rich context before interpreting clustered or duplicated candidates.

## Domain And Annotation-Rescue Rules

- HMMER/Pfam domain support means a sequence contains a recognizable conserved domain under configured thresholds; it does not prove biological function.
- `NO_EXPECTED_DOMAIN_DETECTED` means domain evidence was not detected under current settings; it is not gene absence or loss.
- InterProScan may be used later as broader independent annotation support, but it should not replace locus-level validation for duplication claims.
- Protein-to-genome alignments from miniprot are candidate evidence. They require manual review of coordinates, exon/CDS structure, frameshift/stop flags, and paralog identity before biological interpretation.
- Do not call pseudogenization, loss of function, or gene-state changes from an automated alignment alone.
- For `TP53`, treat p53-family paralogy and divergent/partial alignments as high-risk until additional queries, domain checks, and cross-resource validation are complete.

## RNA-seq And Expression Rules

- Use public RNA-seq only when metadata are clear enough to support the intended question.
- Do not run or interpret differential expression without a defensible design, sample metadata, replicates, and explicit statistical model.
- Use cautious expression language:
  - `detected`,
  - `expressed in tissue X`,
  - `not detected under these conditions`,
  - `consistent with`,
  - `exploratory`.
- Do not use `activated` unless the transcriptomic design supports pathway activation inference.
- Tissue-specific findings, especially retina, must not be generalized to whole-organism longevity without independent support.

## Comparative Longevity Interpretation Rules

- Treat extreme longevity as a multifactorial phenotype. Do not attribute longevity to a single gene, pathway, repeat class, or variant without strong convergent evidence.
- Separate longevity hypotheses from cold/deep-sea adaptation, body size, metabolic rate, population history, ecology, and tissue-specific preservation.
- Do not claim human translational relevance except in a clearly marked speculative section.
- Published claims remain external claims until reproduced or directly supported by this repository's local evidence tables.
- Resource-quality observations such as assembly size, N50, GC content, and BUSCO scores are not biological mechanism evidence.

## Reporting And Figure Rules

- Figures should be generated from versioned result tables, not hand-edited values.
- Diagnostic plots may appear in earlier phases, but report-level figures belong in Phase 9 after tables and schemas are stable.
- Every biological claim in a report should be traceable to `docs/claims_register.md` and to one or more supporting result files.
- Prefer conservative wording in reports:
  - `candidate`,
  - `consistent with`,
  - `plausible lead`,
  - `requires validation`,
  - `artifact-prone`,
  - `not assessed`.
- Avoid promotional language and avoid implying causality, function, adaptation, or clinical relevance without direct evidence.
