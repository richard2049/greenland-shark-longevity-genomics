# Phase 5b Repeat Annotation Setup

Phase 5b is an optional de novo repeat-annotation workflow for the Tokyo/PNAS Greenland shark assembly. It exists to overcome the current limitation that the local Figshare package does not provide parseable repeat annotations.

This workflow is deliberately not part of the default Snakemake run. Repeat annotation is heavy, parameter-sensitive, and best run in Docker or WSL/Linux with explicit provenance.

## Scientific Scope

Allowed:

- Build a species-specific repeat library with RepeatModeler2.
- Annotate repeats with RepeatMasker using the species-specific library.
- Record repeat annotation provenance, checksums, and expected output paths.
- Import the RepeatMasker GFF into Phase 5 and intersect it with Phase 4e candidate loci.
- Treat overlaps as artifact/context evidence.

Not allowed from Phase 5b alone:

- Claims that repeats caused candidate gene duplication.
- Claims of adaptive repeat expansion.
- Claims of longevity mechanism.
- Claims that local repeat overlap proves or disproves a candidate locus.

## Current Preflight

Generate the planning/provenance tables:

```powershell
conda activate green_shark
cd D:\Documents\Greenland_Shark

python -m greenland_shark_longevity.phase5b_repeat_annotation `
  --config config/config.yaml `
  --preflight-output results/repeats/phase5b_repeat_annotation_preflight.tsv `
  --plan-output results/repeats/phase5b_repeat_annotation_plan.tsv `
  --provenance-output data/metadata/phase5b_repeat_annotation_provenance.tsv `
  --output-inventory data/metadata/phase5b_repeat_output_inventory.tsv
```

The current preflight confirms that the Tokyo genome FASTA is available and readable as gzip. RepeatModeler and RepeatMasker outputs are expected but not present until the external run is completed.

## Docker/WSL Strategy

Use Docker Desktop with the `desktop-linux` context, WSL/Linux, or another Linux host. Do not run RepeatModeler/RepeatMasker through native Windows Snakemake.

The current Docker route uses `dfam/tetools:latest` with the resolved digest and observed tool versions recorded in `config/config.yaml` and `data/metadata/phase5b_repeat_annotation_provenance.tsv`. In this image, the relevant executables are not on `PATH`; use:

- `/opt/RepeatModeler/BuildDatabase`
- `/opt/RepeatModeler/RepeatModeler`
- `/opt/RepeatMasker/RepeatMasker`

RepeatModeler 2.0.8 uses `-threads`; do not use the deprecated `-pa` option for RepeatModeler. RepeatMasker 4.2.3 still uses `-pa` for parallelism. Before publication-grade interpretation, prefer pinning the container by digest rather than relying on the moving `latest` tag.

The command templates are written to:

- `results/repeats/phase5b_repeat_annotation_plan.tsv`

Use those commands as the source of truth. They stage an uncompressed genome copy under `data/interim/repeats/phase5b/`, build the RepeatModeler database there, copy the RepeatModeler family library to `results/repeats/phase5b/smic_tokyo-families.fa`, and copy RepeatMasker outputs to stable paths.

For long runs, start RepeatModeler as a detached Docker container and write stdout/stderr to `logs/phase5b_repeatmodeler.*.log`. This avoids frozen PowerShell sessions while preserving inspectable progress and exit status.

Expected final files:

- `results/repeats/phase5b/smic_tokyo-families.fa`
- `results/repeats/phase5b/smic_tokyo.repeatmasker.out.gff`
- `results/repeats/phase5b/smic_tokyo.repeatmasker.out`
- `results/repeats/phase5b/smic_tokyo.repeatmasker.tbl` if produced

Current local status: RepeatModeler produced the species-specific family FASTA, seed alignment, and run log. RepeatMasker produced parseable `.gff` and `.out` files, but no `.tbl`; its log also records a `ProcessRepeats` child process terminated by signal 9. Treat this as sufficient for bounded candidate-locus context after QC, but not as a clean genome-wide repeat summary.

## Import Into Phase 5

After `results/repeats/phase5b/smic_tokyo.repeatmasker.out.gff` exists, add it to `phase5_repeat_context.repeat_annotation_candidates` in `config/config.yaml`. The current config also records the Tokyo assembly report so Phase 5 can map RepeatMasker accessions such as `JBLTJD010000033.1` to Phase 4e names such as `scaffold_33`, and it enables candidate-window filtering:

```yaml
phase5_repeat_context:
  assembly_report: data/raw/references/SMIC_TOKYO_GENOME_2025/GCA_056099535.1_ASM5609953v1_assembly_report.txt
  filter_to_candidate_windows: true
  repeat_annotation_candidates:
    - results/repeats/phase5b/smic_tokyo.repeatmasker.out.gff
```

Then rerun the bounded Phase 5 context step:

```powershell
python -m greenland_shark_longevity.phase5_repeat_context `
  --config config/config.yaml `
  --candidate-loci results/rescue/phase4e_locus_manual_review.tsv `
  --figshare-inventory data/metadata/figshare_annotation_inventory.tsv `
  --resource-status-output data/metadata/phase5_repeat_resource_status.tsv `
  --repeat-features-output results/repeats/phase5_repeat_features.tsv `
  --locus-context-output results/repeats/phase5_candidate_locus_repeat_context.tsv `
  --gene-summary-output results/repeats/phase5_gene_repeat_context_summary.tsv
```

This reuses the existing interval-intersection code and reports repeat context for the `FTH1B`, `H1F0`, `RAD51`, and `TP53` candidate loci. The output `results/repeats/phase5_repeat_features.tsv` is intentionally candidate-window filtered; it is not a genome-wide repeat table.

## Interpretation

Phase 5b can strengthen the repository by replacing `NOT_ASSESSED` repeat context with a reproducible repeat-overlap table. It still cannot make repeat-driven longevity claims. Any overlap should be reported as local context and artifact risk until locus identity, repeat annotation quality, and cross-resource support are evaluated.
