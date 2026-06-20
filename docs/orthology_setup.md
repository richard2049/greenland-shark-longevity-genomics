# Orthology Runtime Setup

Phase 3 uses OrthoFinder for a first-pass, whole-proteome orthogroup screen. The staged inputs are under `data/interim/orthofinder_input/`.

On this native Windows laptop, OrthoFinder is not installed in the `green_shark` conda environment. The Bioconda package is not solvable for the current `win-64` platform because its external search-tool dependency stack is Linux/macOS-oriented. The working strategy is therefore:

1. Run OrthoFinder externally in Docker.
2. Parse the OrthoFinder result directory into repository-standard TSVs.
3. Run `workflow_mode=orthology_postprocess` for candidate mapping and evidence scoring.

## Recommended Route: Docker On Windows

Run from PowerShell in the repository root:

```powershell
$repo = (Get-Location).Path
$run = "results/orthofinder_og_$(Get-Date -Format yyyyMMdd_HHmmss)"

docker --context desktop-linux pull staphb/orthofinder:2.5.5

docker --context desktop-linux run --rm `
  -v "${repo}:/work" `
  -w /work `
  staphb/orthofinder:2.5.5 `
  orthofinder -f data/interim/orthofinder_input -o $run -S diamond -og -t 8 -a 8
```

The `-og` option stops after orthogroup inference. This is enough for the current MVP because the immediate task is candidate copy-number screening, not gene-tree/species-tree inference. DIAMOND is used because it is the practical fast search backend for this small multi-proteome screen.

After OrthoFinder finishes:

```powershell
python -m greenland_shark_longevity.orthofinder_workflow parse `
  --results-dir $run `
  --manifest results/orthology/orthofinder_input_manifest.tsv `
  --gene-count-output results/orthology/orthogroup_gene_counts_long.tsv `
  --species-summary-output results/orthology/orthofinder_species_summary.tsv

snakemake --snakefile workflow/Snakefile --cores 1 --config workflow_mode=orthology_postprocess
```

Expected parsed/postprocessed outputs:

- `results/orthology/orthogroup_gene_counts_long.tsv`
- `results/orthology/orthofinder_species_summary.tsv`
- `results/orthology/candidate_copy_number.tsv`
- `results/orthology/candidate_duplication_audit.tsv`
- `results/evidence/integrated_evidence.tsv`
- `results/validation/candidate_isoform_audit.tsv`
- `results/validation/domain_check_targets.tsv`
- `results/validation/high_priority_rescue_targets.tsv`
- `data/interim/candidate_validation/representative_candidate_proteins.faa`

Record the Docker image tag, command, result directory, thread count, and date in notes or a run log when using this route. The repository captures parsed outputs, but the container invocation itself is external to the Snakemake `orthofinder_run` rule.

## Alternative Route: WSL Or Linux

Run these commands from a WSL/Linux shell, not from native PowerShell:

```bash
cd /mnt/d/Documents/Greenland_Shark
conda env create -f environment.orthology.yml
conda activate green_shark_orthology
orthofinder -h
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=orthology
```

This route lets Snakemake run OrthoFinder directly through the `orthology` workflow mode. Use it when `orthofinder` is available on `PATH`.

## Interpretation

OrthoFinder outputs are genome-wide computational evidence inputs, not final biological claims. Candidate postprocessing maps exact annotation symbols/compact aliases to orthogroups and reports locus-collapsed candidate copy number separately from raw OrthoFinder protein counts. This is intentional because the focal protein FASTA contains isoforms.

The validation-preparation step performs isoform filtering by selecting the longest protein isoform per mapped locus and prepares representative sequences for external domain tools. Domain integrity, cross-resource Greenland shark replication, and candidate-specific rescue are still required before duplication or mechanism-level interpretation.

## Domain Validation After Orthology

After `workflow_mode=orthology_postprocess` has produced `data/interim/candidate_validation/representative_candidate_proteins.faa`, run HMMER/Pfam as the first domain-validation layer:

```bash
hmmpress data/raw/references/PFAM/Pfam-A.hmm
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=domain_validation
```

The domain workflow writes `results/validation/hmmer_pfam_preflight.tsv`, raw HMMER output under `results/domains/`, and parsed classifications in `results/validation/candidate_domain_integrity.tsv`.

The parser uses HMMER `--domtblout` because it records protein identifiers, Pfam accessions, independent domain e-values, bit scores, HMM coordinates, protein coordinates, and domain coverage in a stable machine-readable format. The current classification is intentionally conservative:

- `DOMAIN_SUPPORTED` means at least one Pfam hit passes the configured e-value and full-domain HMM coverage thresholds.
- `PARTIAL_DOMAIN` means a significant hit is present but coverage is incomplete.
- `NO_EXPECTED_DOMAIN_DETECTED` means no qualifying Pfam hit was detected under the configured thresholds.
- `NOT_ASSESSED` means the protein or HMMER/Pfam run is unavailable.

These labels are domain-evidence labels only. They do not prove functional activity, inactivation, loss of function, true duplication, or longevity relevance. InterProScan can be added later as a broader independent annotation layer after this narrower HMMER/Pfam check is stable.

## Annotation Rescue After Domain Validation

After `workflow_mode=domain_validation`, run the targeted Phase 4 rescue layer for unresolved high-priority genes:

```bash
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=annotation_rescue
```

This mode uses reciprocal `phmmer` protein similarity to prioritize candidate focal proteins for genes that exact-symbol mapping left unresolved. It is intentionally limited to `H1F0`, `FTH1B`, `TP53`, and `RAD51` in the current configuration. The output is a protein-level rescue table and a genome-rescue plan; it is not a validated annotation or duplication result.

The genome-level gate is implemented as:

```bash
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=genome_validation
```

This mode prepares miniprot queries from Phase 4 results and runs miniprot only when the configured genome FASTA and executable are locally available. On the current Windows laptop, miniprot is provided by `tools/miniprot-docker.cmd`, and the configured profile uses `candidate_scaffolds` to avoid full-genome Docker memory failure. The resulting alignments are candidate genome-validation inputs and still require manual exon/locus/domain review.
