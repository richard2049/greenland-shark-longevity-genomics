# Greenland Shark Longevity Genomics

This repository is an MVP workflow for a conservative reproducible analysis of public genomic and transcriptomic evidence related to extreme longevity and long-term tissue maintenance in the Greenland shark (*Somniosus microcephalus*).

The workflow is designed to rank evidence, not to overclaim mechanisms. Public-resource claims from papers are registered as external context until reproduced by this repository.

## What This MVP Does

- Registers verified public Greenland shark genome and retinal RNA-seq resources before any large download.
- Validates resource and species metadata into machine-readable TSV files.
- Validates curated candidate gene panels across priority mechanisms.
- Runs assembly/protein/annotation QC on demo fixtures and selected public reference files.
- Stages a small chondrichthyan OrthoFinder input set from the focal Greenland shark protein FASTA and four comparator proteomes.
- Parses externally generated OrthoFinder orthogroups and maps curated candidates onto exact GFF symbols/compact aliases.
- Produces schema-valid copy-number, duplication-audit, evidence-tier, and preliminary report outputs.
- Runs a targeted Phase 4 protein-level annotation-rescue screen for unresolved high-priority genes marked as annotation uncertainty.
- Runs a limited-scope Phase 4b miniprot genome-validation gate for rescued candidate proteins using candidate scaffolds from the local Tokyo genome FASTA.
- Runs a limited Phase 5 repeat-context check around hardened candidate loci using existing repeat annotations or externally generated RepeatMasker annotations, with NCBI accession-to-scaffold normalization and candidate-window filtering.
- Provides an optional Phase 5b planning/provenance layer for de novo repeat annotation in Docker or WSL/Linux.
- Registers Phase 7 RNA-seq metadata and creates a candidate-expression quantification strategy without downloading or quantifying raw reads.
- Constructs a Phase 7c annotation-derived transcript FASTA and `tx2gene` table, then validates candidate IDs before any raw-read quantification.
- Adds Phase 7d retina-only raw-read intake, lightweight FASTQ QC, and guarded Salmon quantification when paired FASTQs and a Salmon index are supplied locally.
- Adds Phase 7e expression-interpretation hardening so low mapping rates, reference ambiguity, duplicated candidate loci, and repeat-context caveats are reviewed before Phase 8 evidence scoring.
- Adds Phase 8a expression-evidence integration so retina expression support is carried into evidence tables without automatic biological tier upgrades.
- Adds Phase 8b final integrated evidence scoring and tier audit before Phase 9 reporting.

The demo fixtures are marked `DEMO_ONLY_NOT_BIOLOGICAL_EVIDENCE` and must not be interpreted biologically.

## What This MVP Cannot Claim

- It does not infer gene loss from annotation absence.
- It does not claim functional activation, inactivation, causal longevity effects, telomere length, telomerase activity, repeat-driven adaptation, or human translational relevance.
- It does not treat OrthoFinder orthogroups alone as validated gene duplication, protein-domain integrity, or functional evidence.
- It does not run BUSCO, selection tests, InterProScan, or full-genome miniprot rescue yet.
- It does not run optional Phase 5b de novo repeat annotation automatically; RepeatModeler/RepeatMasker must be launched manually in a suitable Linux/container environment.
- It does not treat Phase 4 protein-level rescue candidates or Phase 4b miniprot alignments as final gene annotations or validated duplications.
- It does not treat repeat overlap as evidence for repeat-mediated causation, adaptation, functional effect, or longevity relevance.
- It does not treat Phase 7b quantification planning as evidence that a candidate gene is detected or expressed.
- It does not treat Phase 7c transcript-reference presence as expression evidence.
- It does not download raw RNA-seq reads automatically, run differential expression, infer pathway activity, or generalize retina quantification to whole-organism aging.
- It does not let Phase 7e/8a retina expression support or Phase 8b reporting classes upgrade biological evidence tiers unless upstream validation criteria are satisfied.
- It does not treat published claims as reproduced until local evidence tables support them.

## Verified Public Sources Registered Initially

| Resource | Accession | Usage in MVP |
|---|---:|---|
| FLI/University of Copenhagen Greenland shark genome | `GCA_052327205.1_FLI_Smic_1.0`, `PRJNA1236315`, `SAMN47384478` | Manifest only |
| University of Tokyo Somniosus sequencing BioProject | `PRJNA1218902`, `GCA_056099535.1`, `JBLTJD000000000` | Manifest only |
| Retinal/genomic read study | `PRJNA1246101` | Manifest only |

## Quickstart

Create the environment:

```bash
conda env create -f environment.yml
conda activate green_shark
```

Run the default demo workflow:

```bash
snakemake --snakefile workflow/Snakefile --cores 1
```

On Windows systems where the default user cache directory is locked, run Snakemake with a workspace-local source cache:

```bash
snakemake --snakefile workflow/Snakefile --cores 1 \
  --shared-fs-usage input-output persistence software-deployment software-deployment-cache sources storage-local-copies \
  --runtime-source-cache-path .snakemake/source-cache-runtime
```

Run only metadata/provenance:

```bash
snakemake --snakefile workflow/Snakefile --cores 1 --config workflow_mode=metadata_only
```

Run Phase 1 reference intake and resource-quality QC:

```bash
snakemake --snakefile workflow/Snakefile --cores 1 --config workflow_mode=reference_only
```

This downloads only selected small NCBI files by default: assembly statistics, assembly reports, feature-count files, checksums, and the `PRJNA1246101` coding-sequence/CDS-translation FASTA files. Genome FASTA and GenBank flatfiles are registered in `data/metadata/reference_file_inventory.tsv` but are not downloaded by default because they are gigabyte-scale. Separate protein FASTA and genomic GFF files were not listed in the registered NCBI assembly directories during the 2026-05-27 intake pass, so assembly-package protein sequence QC remains `NOT_ASSESSED`.

The reference workflow also audits registered publication claims in `config/publication_claims.yaml`. PNAS/Yang et al. 2026 claims about linker histone H1.0 and `FTH1B` copy number are tracked as `NOT_ASSESSED` until the repository has genome-wide protein/annotation resources and orthology/domain/locus validation.

The Phase 3 readiness step records the search for PNAS/bioRxiv/NCBI/author-linked protein FASTA, GFF/GTF, gene-model, and orthogroup/copy-number inputs in `results/readiness/phase3_source_checks.tsv`. Manual inspection of the PNAS data availability section identified a Figshare annotation source at `https://figshare.com/s/4f1adabbc84fbf5a72e0`, including protein sequences predicted from transcriptome data. The manually downloaded files are inventoried in `data/metadata/figshare_annotation_inventory.tsv`.

The workflow now treats `complete.proteins.faa` as the candidate genome-annotation protein FASTA for the Tokyo/PNAS resource, runs real protein QC on it, parses coordinates from `complete.genomic.gff`, and stages `data/interim/orthofinder_input/smic__SMIC_TOKYO_GENOME_2025.faa`. This is still input readiness, not biological evidence. The transcriptome-derived TransDecoder peptide FASTA is inventoried but kept separate from genome-wide copy-number analysis. The small `PRJNA1246101` CDS translations remain excluded from OrthoFinder because they are not a genome-wide proteome.

The first comparator intake uses RefSeq protein FASTA files for four chondrichthyans: elephant shark (`Callorhinchus milii`), small-spotted catshark (`Scyliorhinus canicula`), thorny skate (`Amblyraja radiata`), and whale shark (`Rhincodon typus`). These were selected before more distant outgroups because the immediate goal is to make the Greenland shark candidate-gene analysis comparable within cartilaginous fishes. Their protein FASTA files are downloaded, checksummed, QCed, and staged under `data/interim/orthofinder_input/`.

### Phase 3 Orthology On This Laptop

On this native Windows laptop, the practical route is to run OrthoFinder in Docker, then let this repository parse and postprocess the result. This keeps the heavy external binary out of the Windows conda environment while preserving reproducible downstream tables.

From PowerShell in the repository root:

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

The `-og` option is intentional for this MVP: Phase 3 currently needs orthogroups and gene-count tables, not full gene-tree/species-tree inference. DIAMOND is used because it is the appropriate fast first-pass search backend for a small multi-proteome orthogroup screen.

Then parse and postprocess the result:

```powershell
python -m greenland_shark_longevity.orthofinder_workflow parse `
  --results-dir $run `
  --manifest results/orthology/orthofinder_input_manifest.tsv `
  --gene-count-output results/orthology/orthogroup_gene_counts_long.tsv `
  --species-summary-output results/orthology/orthofinder_species_summary.tsv

snakemake --snakefile workflow/Snakefile --cores 1 --config workflow_mode=orthology_postprocess
```

The current postprocessing maps exact GFF gene symbols or compact aliases onto `Orthogroups.tsv`, collapses protein isoforms to unique Greenland shark loci for candidate copy-number reporting, and keeps raw OrthoFinder protein counts as separate context. This avoids treating isoforms as duplications.

The same `orthology_postprocess` mode also performs candidate-specific validation preparation:

- `results/validation/candidate_isoform_audit.tsv` collapses candidate proteins to locus-level representative isoforms.
- `data/interim/candidate_validation/representative_candidate_proteins.faa` stores representative candidate proteins for downstream domain tools.
- `results/validation/domain_check_targets.tsv` records which proteins are ready for InterProScan or HMMER/Pfam domain checks.
- `results/validation/high_priority_rescue_targets.tsv` records unresolved high-priority genes requiring reciprocal similarity search and protein-to-genome rescue.

This step does not assess protein-domain integrity by itself. Domain calls require an explicit domain tool/database run with recorded database versions, e-values, and coverage.

### Candidate Domain Validation

The next validation layer uses HMMER/Pfam before broader InterProScan annotation. HMMER/Pfam is the first-pass choice because the immediate question is whether representative candidate proteins contain qualifying conserved domains with explicit e-values and coverage. InterProScan remains useful later for broader independent annotation, but it is heavier than needed for the current 24 representative proteins.

On this native Windows environment, the workflow uses the `pyhmmer` backend because command-line `hmmscan` is not available through the current `win-64` conda channels. This is still a HMMER/Pfam profile-HMM scan, but it does not require a separate `hmmscan.exe`.

Configure a local Pfam-A HMM database in `config/config.yaml`:

```yaml
domain_validation:
  backend: pyhmmer
  hmmer_executable: hmmscan
  pfam_hmm: data/raw/references/PFAM/Pfam-A.hmm
  pfam_release: "38.1"
```

For command-line `hmmscan`, press the Pfam database once:

```bash
hmmpress data/raw/references/PFAM/Pfam-A.hmm
```

The current `pyhmmer` backend can read the plain `Pfam-A.hmm` file directly; pressed indexes are not required for this backend.

Then run the domain-validation mode:

```bash
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=domain_validation
```

Expected domain-validation outputs:

- `results/validation/hmmer_pfam_preflight.tsv`
- `results/domains/representative_candidate_proteins.pfam.domtblout`
- `results/domains/representative_candidate_proteins.pfam.txt`
- `results/validation/candidate_domain_integrity.tsv`

The final table classifies representative proteins as `DOMAIN_SUPPORTED`, `PARTIAL_DOMAIN`, `NO_EXPECTED_DOMAIN_DETECTED`, or `NOT_ASSESSED`. These are protein-domain validation statuses only. `NO_EXPECTED_DOMAIN_DETECTED` must not be interpreted as gene absence, loss of function, or biological inactivation.

See `docs/domain_validation.md` for setup and interpretation notes.

### Phase 4 Targeted Annotation Rescue

Phase 4 is focused only on unresolved high-priority genes currently marked as `ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH`: `H1F0`, `FTH1B`, `TP53`, and `RAD51`.

The implemented first pass uses `pyhmmer`/HMMER `phmmer` reciprocal protein similarity. This is the appropriate first screen because the repository already has staged comparator proteomes and a focal protein FASTA. Product-description pattern filters are intentionally gene-specific and conservative to avoid merging nearby paralogs such as H1 variants, RAD51 paralogs, or ferritin-family members.

Run:

```bash
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=annotation_rescue
```

Expected Phase 4 outputs:

- `results/rescue/phase4_rescue_query_inventory.tsv`
- `results/rescue/phase4_forward_hits.tsv`
- `results/rescue/phase4_reciprocal_hits.tsv`
- `results/rescue/phase4_annotation_rescue_summary.tsv`
- `results/rescue/phase4_genome_rescue_plan.tsv`
- `data/interim/annotation_rescue/phase4_rescue_candidate_proteins.faa`

Current Phase 4 interpretation is deliberately narrow: reciprocal protein-level rescue identified candidate focal protein models for `H1F0`, `FTH1B`, and `RAD51`, while `TP53` did not produce a reciprocal protein-level rescue candidate under the configured filters. This is not evidence that `TP53` is absent. All four genes still require genome-level rescue or locus inspection before annotation, duplication, loss, or mechanism-level language is justified.

See `docs/annotation_rescue.md` for method and interpretation details.

Run the Phase 4b genome-validation gate:

```bash
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=genome_validation
```

This prepares miniprot query FASTA/metadata and runs protein-to-genome alignment only when the configured genome FASTA and `miniprot` executable are available. On this Windows laptop, `miniprot` is provided through `tools/miniprot-docker.cmd`, and `config/config.yaml` uses the local Tokyo genome FASTA in `candidate_scaffolds` mode. This mode extracts only scaffolds already implicated by Phase 4 protein-level rescue candidates, avoiding a full-genome Docker memory failure while preserving a conservative validation gate.

Current Phase 4b output should be read narrowly. The workflow detects separable high-coverage candidate-scaffold alignments for `H1F0`, `FTH1B`, and `RAD51`; `TP53` remains `not_assessable` because it has no protein-level rescued candidate region under the limited scaffold scope. These are genome-alignment validation inputs, not evidence of function, adaptation, loss, inactivation, or validated duplication.

Run the Phase 4c/4d review and consolidation steps after Phase 4b outputs are current:

```bash
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=phase4d_consolidation
```

Phase 4c reviews miniprot loci, Pfam support for rescue queries, overlap among candidate loci, and a bounded targeted `TP53` search against scaffolds selected from Phase 4 forward-hit coordinates. Phase 4d then writes `results/evidence/phase4d_candidate_interpretation.tsv` and updates `results/evidence/integrated_evidence.tsv` while preserving the pre-consolidation score in `results/evidence/phase3_integrated_evidence.tsv`.

Current Phase 4d interpretation is conservative: `H1F0`, `FTH1B`, and `RAD51` are Tier 2 candidate-locus findings that require manual locus review and cross-resource support before any duplication or function language is justified. `TP53` remains `Artifact/uncertain`; the targeted p53-family alignment is not a gene-state or mechanism claim.

Run the Phase 4e manual-review hardening layer when Phase 4c/4d outputs and the local Figshare annotation GFF are present. From PowerShell:

```powershell
python -m greenland_shark_longevity.phase4e_locus_hardening `
  --phase4c-locus-review results/rescue/phase4c_locus_review.tsv `
  --domain-integrity results/rescue/phase4c_rescue_domain_integrity.tsv `
  --tp53-alignment-hits results/rescue/tp53_targeted_forward_alignment_hits.tsv `
  --tp53-target-regions results/rescue/tp53_forward_hit_target_regions.tsv `
  --annotation-gff "data/raw/references/SMIC_TOKYO_GENOME_2025/figshare_annotation/greenland shark annotation/complete.genomic.gff" `
  --locus-output results/rescue/phase4e_locus_manual_review.tsv `
  --summary-output results/rescue/phase4e_gene_hardened_summary.tsv `
  --evidence-output results/evidence/phase4e_hardened_evidence.tsv
```

Phase 4e records exact/focal annotation overlap, product-consistent family overlap, candidate-locus spacing, domain support, and artifact flags. Current Phase 4e output keeps `H1F0` and `RAD51` as Tier 2 candidate loci requiring cross-resource validation, keeps `FTH1B` as a high-risk ferritin-family cluster without validated copy-number interpretation, and keeps `TP53` as `Artifact/uncertain`.

### Phase 5 Limited Repeat Context

Phase 5 is implemented only in candidate-locus mode. It does not run RepeatMasker, RepeatModeler, or de novo repeat discovery inside the default workflow. Instead, it checks the current Tokyo/PNAS Figshare inventory and any explicitly configured repeat-annotation paths, normalizes RepeatMasker GenBank accessions to assembly-report scaffold names, filters repeat records to configured candidate windows, and intersects them with the Phase 4e loci for:

- `FTH1B` ferritin-family cluster on `scaffold_33`
- `H1F0` candidate loci on `scaffold_33`
- `RAD51` candidate loci on `scaffold_12`
- `TP53` targeted p53-family locus on `scaffold_247`

Current Phase 5 output imports the externally generated `results/repeats/phase5b/smic_tokyo.repeatmasker.out.gff` and retains 1,204 repeat features near the Phase 4e candidate windows. The full RepeatMasker run produced parseable `.gff` and `.out` files but no `.tbl`, and its log contains a `ProcessRepeats` signal-9 warning. Therefore these outputs are appropriate for candidate-locus artifact/context review, not genome-wide repeat-summary claims.

Run the bounded Phase 5 step from PowerShell:

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

After running the Phase 5c QC step below, integrate Phase 5 into the evidence/reporting layer as artifact risk only:

```powershell
python -m greenland_shark_longevity.phase5_evidence_integration `
  --base-evidence results/evidence/phase4d_integrated_evidence.tsv `
  --phase4e-evidence results/evidence/phase4e_hardened_evidence.tsv `
  --phase5-gene-summary results/repeats/phase5_gene_repeat_context_summary.tsv `
  --phase5c-gene-qc results/repeats/phase5c_gene_repeat_qc_summary.tsv `
  --phase5-evidence-output results/evidence/phase5_repeat_context_evidence.tsv `
  --integrated-output results/evidence/integrated_evidence.tsv `
  --report-output reports/generated/phase5_repeat_context_report.md
```

This step does not upgrade evidence tiers. It adds local repeat-context artifact-risk caveats to `FTH1B`, `H1F0`, `RAD51`, and `TP53`, with direct repeat-overlap risk especially relevant for `FTH1B`, `RAD51`, and `TP53`.

### Phase 5b Optional Repeat Annotation

Phase 5b is the path for overcoming the current missing-repeat-annotation limitation. It prepares a reproducible external RepeatModeler2/RepeatMasker route, but does not run the heavy tools from the default workflow.

Generate the Phase 5b planning/provenance tables:

```powershell
python -m greenland_shark_longevity.phase5b_repeat_annotation `
  --config config/config.yaml `
  --preflight-output results/repeats/phase5b_repeat_annotation_preflight.tsv `
  --plan-output results/repeats/phase5b_repeat_annotation_plan.tsv `
  --provenance-output data/metadata/phase5b_repeat_annotation_provenance.tsv `
  --output-inventory data/metadata/phase5b_repeat_output_inventory.tsv
```

Expected Phase 5b planning outputs:

- `results/repeats/phase5b_repeat_annotation_preflight.tsv`
- `results/repeats/phase5b_repeat_annotation_plan.tsv`
- `data/metadata/phase5b_repeat_annotation_provenance.tsv`
- `data/metadata/phase5b_repeat_output_inventory.tsv`

The recommended external route is Docker `desktop-linux` or WSL/Linux. The current Docker route uses `dfam/tetools:latest` with its resolved digest, RepeatModeler 2.0.8, and RepeatMasker 4.2.3 recorded in config/provenance. In that image, use `/opt/RepeatModeler/BuildDatabase`, `/opt/RepeatModeler/RepeatModeler`, and `/opt/RepeatMasker/RepeatMasker`; RepeatModeler uses `-threads`, while RepeatMasker uses `-pa`. See `docs/repeat_annotation_setup.md`.

### Phase 5c Repeat-Context QC Hardening

Phase 5c hardens the candidate-locus repeat context without rerunning RepeatModeler or RepeatMasker. It streams the large RepeatMasker `.out` file, keeps only records intersecting the Phase 4e candidate windows, compares `.out` direct-overlap calls against the GFF-derived Phase 5 table, records run-integrity warnings, and computes local window repeat density.

Run it from PowerShell:

```powershell
python -m greenland_shark_longevity.phase5c_repeat_qc `
  --config config/config.yaml `
  --candidate-loci results/rescue/phase4e_locus_manual_review.tsv `
  --phase5-locus-context results/repeats/phase5_candidate_locus_repeat_context.tsv `
  --phase5b-inventory data/metadata/phase5b_repeat_output_inventory.tsv `
  --integrity-output results/repeats/phase5c_repeatmasker_integrity.tsv `
  --locus-qc-output results/repeats/phase5c_locus_repeat_qc.tsv `
  --gene-qc-output results/repeats/phase5c_gene_repeat_qc_summary.tsv `
  --report-output reports/generated/phase5c_repeat_context_qc_report.md
```

Current Phase 5c status:

- RepeatMasker `.gff` and `.out` are present and parseable.
- RepeatMasker `.tbl` is missing, so genome-wide repeat percentages are not supported.
- Logs retain a `ProcessRepeats` signal-9 warning and an LTRPipeline warning, so outputs remain `candidate_locus_context_only`.
- `.out` parsing supports the direct-overlap calls for `FTH1B`, `RAD51`, and `TP53`; `H1F0` remains window-context only.

### Phase 6 Telomere Motif And Gene Readiness

Phase 6 is implemented as a conservative telomere-readiness layer. It searches exact canonical vertebrate telomeric motifs (`TTAGGG` and `CCCTAA`) in scaffold-end windows and audits curated telomere/shelterin genes against the current candidate/evidence tables.

Run it from PowerShell:

```powershell
python -m greenland_shark_longevity.phase6_telomere `
  --config config/config.yaml `
  --motif-scan-output results/telomere/phase6_telomeric_motif_scan.tsv `
  --enrichment-output results/telomere/phase6_scaffold_end_enrichment.tsv `
  --gene-audit-output results/telomere/phase6_telomere_gene_audit.tsv `
  --report-output reports/generated/phase6_telomere_report.md
```

Current Phase 6 status:

- 2,136 Tokyo assembly sequences were scanned with 10 kb terminal windows.
- Terminal motif density is higher than adjacent internal control windows in this exact-motif scan (`end_enrichment_ratio = 1.61675`).
- `TERT`, `TERF1`, `TERF2`, `POT1`, `TINF2`, `ACD`, and `TERF2IP` currently map as single-copy exact-symbol candidates in the post-OrthoFinder candidate table.
- These are sequence/context and readiness observations only. They do not support telomere length, telomerase activity, rejuvenation, pathway activity, or longevity-mechanism claims.

### Phase 7a RNA-seq Metadata And Expression Readiness

Phase 7a is implemented as a metadata/readiness layer for `PRJNA1246101`. It records exact SRA run accessions and classifies which runs are suitable for a future candidate-expression audit. It does not download raw reads, quantify expression, test differential expression, or report expression support.

Run it from PowerShell:

```powershell
python -m greenland_shark_longevity.phase7_rnaseq `
  --config config/config.yaml `
  --manifest-output data/metadata/rnaseq_manifest.tsv `
  --readiness-output results/rnaseq/phase7_rnaseq_readiness.tsv `
  --expression-plan-output results/rnaseq/phase7_candidate_expression_plan.tsv `
  --report-output reports/generated/phase7_rnaseq_readiness_report.md
```

Current Phase 7a status:

- Four `PRJNA1246101` SRA runs are registered from NCBI RunInfo.
- Three PolyA paired-end retinal RNA-seq runs are metadata-ready for future candidate expression audit: `SRR32965275`, `SRR32965277`, and `SRR32965276`.
- One paired-end WGS run, `SRR32965274`, is retained for provenance but excluded from expression planning.
- All 41 curated candidate-panel genes are staged for future retina-only expression audit with `NOT_QUANTIFIED_PHASE7A_METADATA_ONLY` status.

### Phase 7b Candidate Expression Quantification Strategy

Phase 7b defines the future candidate-expression quantification route and maps candidate-panel genes onto available annotation/CDS identifiers. It does not download reads, run Salmon, quantify expression, test differential expression, or report expression support.

The recommended first-pass method is Salmon-style transcript quantification from a Tokyo annotation-derived transcript FASTA with genome decoys and a transcript-to-gene map. Salmon is the preferred first pass because the current expression task is limited to three paired-end retinal RNA-seq runs and candidate-level detection/expression summaries; it is faster and easier to audit than a full genome-alignment workflow for every run. Genome-aligned counting with STAR or HISAT2 plus featureCounts remains the validation route for ambiguous duplicated/paralogous loci such as `FTH1B`, `H1F0`, `RAD51`, and `TP53`.

Run the strategy layer from PowerShell:

```powershell
python -m greenland_shark_longevity.phase7b_quantification_strategy `
  --config config/config.yaml `
  --reference-strategy-output results/rnaseq/phase7b_reference_quantification_strategy.tsv `
  --candidate-quant-map-output results/rnaseq/phase7b_candidate_quantification_map.tsv `
  --run-plan-output results/rnaseq/phase7b_quantification_run_plan.tsv `
  --report-output reports/generated/phase7b_quantification_strategy_report.md
```

Current Phase 7b status:

- All 41 curated candidate-panel genes are evaluated against the current Phase 7a plan, candidate copy-number table, isoform audit, Phase 4e hardening tables, and local `complete.cds.fna` headers.
- The output distinguishes candidates ready for first-pass retina quantification from candidates that require targeted locus or paralog review before expression-support interpretation.
- `FTH1B`, `H1F0`, `RAD51`, and `TP53` remain high-ambiguity candidates; future expression summaries for these genes must include multi-mapping/locus-review caveats.
- Phase 7b outputs are quantification planning and candidate-reference mapping only. They do not show that any gene is detected or expressed.

### Phase 7c Expression Reference Construction

Phase 7c constructs the annotation-derived transcript FASTA and `tx2gene` table needed for future Salmon quantification. It also validates whether Phase 7b candidate gene IDs are represented in that reference. It does not download raw reads, build a Salmon index, quantify expression, or report expression support.

The implemented extractor is pure Python and streams the genome FASTA one scaffold at a time. This is intentional: the Tokyo assembly is multi-gigabase, and loading the whole genome into memory would be fragile on this laptop. The script uses GFF3 `exon` intervals for transcript reconstruction and records any CDS-only fallback separately in the QC table.

Run it from PowerShell:

```powershell
python -m greenland_shark_longevity.phase7c_expression_reference `
  --config config/config.yaml `
  --transcript-fasta-output data/interim/rnaseq/phase7c/smic_tokyo.annotation_transcripts.fna `
  --tx2gene-output data/interim/rnaseq/phase7c/smic_tokyo.tx2gene.tsv `
  --reference-qc-output results/rnaseq/phase7c_expression_reference_qc.tsv `
  --candidate-validation-output results/rnaseq/phase7c_candidate_reference_validation.tsv `
  --report-output reports/generated/phase7c_expression_reference_report.md
```

Expected Phase 7c outputs:

- `data/interim/rnaseq/phase7c/smic_tokyo.annotation_transcripts.fna`
- `data/interim/rnaseq/phase7c/smic_tokyo.tx2gene.tsv`
- `results/rnaseq/phase7c_expression_reference_qc.tsv`
- `results/rnaseq/phase7c_candidate_reference_validation.tsv`
- `reports/generated/phase7c_expression_reference_report.md`

Phase 7c is still an input-validation layer. A transcript being present in the reference does not mean it is expressed in retina.

### Phase 7d Retina RNA-seq Intake/QC And Salmon Quantification

Phase 7d is the first execution-capable RNA-seq quantification layer. It checks for local paired retinal FASTQs, records file sizes and missing-read status, performs lightweight FASTQ sanity QC when files are present, checks Salmon/index readiness, and runs Salmon only when configured inputs are available. It does not download SRA reads automatically, run differential expression, infer pathway activity, or claim organism-wide aging relevance.

Salmon is used as the first-pass quantifier because the current expression objective is candidate-level retina support from a small paired-end dataset. It is efficient, records transcript-level ambiguity, and can be summarized through `tx2gene`. Genome-aligned counting remains the follow-up validation path for duplicated or paralogous candidates such as `FTH1B`, `H1F0`, `RAD51`, and `TP53`.

On this native-Windows laptop, the configured Salmon executable is `tools/salmon-docker.cmd`, which runs `combinelab/salmon:latest` through Docker Desktop and translates repository-relative paths for the Linux container. This is the preferred local route because Salmon is Linux-first in many bioinformatics environments and the repo already uses containerized execution for tools that are unreliable on native Windows.

Run it from PowerShell after local FASTQs and, if quantification is intended, a Salmon index are available:

```powershell
python -m greenland_shark_longevity.phase7d_rnaseq_quantification `
  --config config/config.yaml `
  --raw-read-intake-output results/rnaseq/phase7d_raw_read_intake.tsv `
  --fastq-qc-output results/rnaseq/phase7d_fastq_qc.tsv `
  --salmon-preflight-output results/rnaseq/phase7d_salmon_preflight.tsv `
  --salmon-quant-summary-output results/rnaseq/phase7d_salmon_quant_summary.tsv `
  --candidate-expression-matrix-output results/rnaseq/phase7d_candidate_expression_matrix.tsv `
  --report-output reports/generated/phase7d_rnaseq_quantification_report.md
```

Expected Phase 7d outputs:

- `results/rnaseq/phase7d_raw_read_intake.tsv`
- `results/rnaseq/phase7d_fastq_qc.tsv`
- `results/rnaseq/phase7d_salmon_preflight.tsv`
- `results/rnaseq/phase7d_salmon_quant_summary.tsv`
- `results/rnaseq/phase7d_candidate_expression_matrix.tsv`
- `reports/generated/phase7d_rnaseq_quantification_report.md`

If reads or Salmon are not available, Phase 7d writes explicit `NOT_RUN` status rows rather than silently skipping or inventing expression support. If a previous run already produced a non-empty `quant.sf`, Phase 7d reuses that file as a restart checkpoint; delete that sample's `results/rnaseq/salmon/<run>/` directory only when intentional re-quantification is required.

### Phase 7e Expression Interpretation Hardening

Phase 7e reviews the Phase 7d expression matrix before expression support can be used in Phase 8. It joins Salmon mapping metrics, lightweight FASTQ QC, Phase 7c reference ambiguity, Phase 3/4 duplication and locus-review outputs, and Phase 5 repeat-context artifact flags. It does not rerun quantification, test differential expression, infer pathway activity, or claim whole-organism aging relevance.

The implemented method is deterministic table integration rather than a new statistical model. That is intentional: the current dataset has three retina RNA-seq runs and no defensible contrast for differential expression. Salmon remains the right first-pass quantifier for this stage, but ambiguous loci such as `FTH1B`, `H1F0`, `RAD51`, and `TP53` are deferred until locus/reference ambiguity is resolved.

Run it from PowerShell after Phase 7d has completed:

```powershell
python -m greenland_shark_longevity.phase7e_expression_hardening `
  --config config/config.yaml `
  --run-qc-output results/rnaseq/phase7e_run_qc_review.tsv `
  --candidate-output results/rnaseq/phase7e_candidate_expression_hardened.tsv `
  --parameter-review-output results/rnaseq/phase7e_parameter_review.tsv `
  --report-output reports/generated/phase7e_expression_hardening_report.md
```

Expected Phase 7e outputs:

- `results/rnaseq/phase7e_run_qc_review.tsv`
- `results/rnaseq/phase7e_candidate_expression_hardened.tsv`
- `results/rnaseq/phase7e_parameter_review.tsv`
- `reports/generated/phase7e_expression_hardening_report.md`

Current parameter policy: `TPM >= 1` plus `NumReads >= 10` is treated as an exploratory detected/not-detected threshold, not a statistical model. At least two of three retina runs are required for consistent low-ambiguity tissue support. Mapping-rate thresholds flag review needs only: `<30%` is cautious use and `<20%` requires review. These thresholds do not discard data or create biological claims.

### Phase 8a Expression Evidence Integration

Phase 8a integrates Phase 7e hardened expression support into the evidence-scoring layer. The method is deterministic rule-based table integration, not a new statistical model. This is the appropriate choice here because the RNA-seq evidence is retina-only, the design does not support differential expression, and several candidates have reference, paralog, or locus ambiguity.

Phase 8a can append retina-specific expression support or expression caveats to candidate evidence rows. It does not upgrade the integrated biological evidence tier. A stronger tier requires independent validation such as orthology support, domain integrity, separable loci where duplication is claimed, cross-resource support, and expression-reference/locus review.

Run it after Phase 7e:

```powershell
python -m greenland_shark_longevity.phase8a_expression_evidence_integration `
  --config config/config.yaml `
  --expression-evidence-output results/evidence/phase8a_expression_support_evidence.tsv `
  --integration-audit-output results/evidence/phase8a_expression_integration_audit.tsv `
  --integrated-output results/evidence/phase8a_integrated_evidence.tsv `
  --report-output reports/generated/phase8a_expression_evidence_report.md
```

Expected Phase 8a outputs:

- `results/evidence/phase8a_expression_support_evidence.tsv`
- `results/evidence/phase8a_expression_integration_audit.tsv`
- `results/evidence/phase8a_integrated_evidence.tsv`
- `reports/generated/phase8a_expression_evidence_report.md`

### Phase 8b Final Evidence Scoring And Tier Audit

Phase 8b consolidates the current evidence streams into the final pre-reporting evidence table. It audits orthology/candidate mapping, Phase 4e locus hardening, Phase 5c repeat-context QC, Phase 6 telomere-gene readiness, and Phase 7e/8a retina expression support.

The implemented method is deterministic rule-based auditing rather than a numeric weighted score. This is intentional: the evidence streams are heterogeneous and incomplete, and a numeric score would imply unsupported precision. Phase 8b may retain or conservatively downgrade tiers, but it does not upgrade tiers. Candidate rows are assigned reporting classes such as `PLAUSIBLE_LEAD_REQUIRES_VALIDATION`, `PLAUSIBLE_LEAD_ARTIFACT_PRONE`, or `ARTIFACT_UNCERTAIN_DO_NOT_CLAIM_BIOLOGICALLY`.

Run it after Phase 8a:

```powershell
python -m greenland_shark_longevity.phase8b_final_evidence_scoring `
  --config config/config.yaml `
  --final-evidence-output results/evidence/phase8b_final_integrated_evidence.tsv `
  --tier-audit-output results/evidence/phase8b_tier_audit.tsv `
  --mechanism-summary-output results/evidence/phase8b_mechanism_summary.tsv `
  --report-output reports/generated/phase8b_final_evidence_report.md
```

Expected Phase 8b outputs:

- `results/evidence/phase8b_final_integrated_evidence.tsv`
- `results/evidence/phase8b_tier_audit.tsv`
- `results/evidence/phase8b_mechanism_summary.tsv`
- `reports/generated/phase8b_final_evidence_report.md`

### Phase 9 Report And Figure Generation

Phase 9 turns the Phase 8b audit tables into a final conservative report and figure set. It uses `phase8b_tier_audit.tsv` and `phase8b_mechanism_summary.tsv` as the primary inputs, with `phase8b_final_integrated_evidence.tsv` retained as the row-level evidence table.

The method is deterministic report generation, not statistical modeling. This is intentional: Phase 8b outputs are categorical evidence-audit classes, not replicate-level measurements with an inferential design. Figures are TSV-backed SVGs with explicit count axes, source captions, interpretation guardrails, and a colorblind-safe palette, so plotted values remain inspectable and reproducible without extra plotting dependencies.

Run it after Phase 8b:

```powershell
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
```

Expected Phase 9 outputs:

- `reports/final/greenland_shark_longevity_phase9_report.md`
- `results/reporting/phase9_figure_manifest.tsv`
- `results/reporting/phase9_key_findings.tsv`
- `results/reporting/phase9_report_package_audit.tsv`
- `results/reporting/phase9_public_repository_readiness.tsv`
- `reports/figures/data/phase9_evidence_tier_summary.tsv`
- `reports/figures/data/phase9_reporting_class_summary.tsv`
- `reports/figures/data/phase9_mechanism_evidence_matrix.tsv`
- `reports/figures/data/phase9_artifact_context_summary.tsv`
- `reports/figures/phase9_evidence_tier_summary.svg`
- `reports/figures/phase9_reporting_class_summary.svg`
- `reports/figures/phase9_mechanism_evidence_matrix.svg`
- `reports/figures/phase9_artifact_context_summary.svg`

### Phase 9 Report Package Audit

After regenerating Phase 9, run the final package audit before public-release planning:

```powershell
python -m greenland_shark_longevity.phase9_package_audit `
  --config config/config.yaml `
  --package-audit-output results/reporting/phase9_report_package_audit.tsv `
  --release-readiness-output results/reporting/phase9_public_repository_readiness.tsv `
  --report-output reports/final/phase9_report_package_audit.md
```

The audit checks report traceability, figure metadata/provenance, language guardrails, and public-repository readiness. It also records a figure-stack decision: TSV-backed SVGs are appropriate for the current categorical evidence-audit figures, while an optional matplotlib/R export layer can be added later for manuscript styling.

The public-release policy is to keep broad `results/` outputs ignored while unignoring only the small final evidence/reporting TSV release bundle by exact `.gitignore` exception. Raw data, intermediates, Salmon outputs, repeat outputs, logs, and large workflow products remain ignored. Use `docs/release_checklist.md` before public release.

Expected audit outputs:

- `results/reporting/phase9_report_package_audit.tsv`
- `results/reporting/phase9_public_repository_readiness.tsv`
- `reports/final/phase9_report_package_audit.md`

### Periodic Output Refresh

On this native Windows laptop, prefer direct Python module commands for periodic refreshes. Broad Snakemake dry-runs have repeatedly hung in this repository, so they are not the routine validation route here.

```powershell
conda activate green_shark
cd D:\Documents\Greenland_Shark

python -m pytest -q `
  tests/test_phase9_package_audit.py `
  tests/test_phase8b_final_evidence_scoring.py `
  tests/test_phase9_report_generation.py `
  tests/test_phase8a_expression_evidence_integration.py `
  tests/test_phase7d_rnaseq_quantification.py `
  tests/test_phase7e_expression_hardening.py `
  tests/test_phase7c_expression_reference.py `
  tests/test_phase7b_quantification_strategy.py `
  tests/test_phase7_rnaseq.py `
  --basetemp .tmp/pytest_phase8_refresh
```

The full direct refresh command sequence through Phase 9 and the package audit is maintained in `docs/runbook.md`.

The closest single Snakemake target for a managed refresh through Phase 9 is available for WSL/Linux, Docker-backed environments, or controlled native-Windows execution when direct modules are insufficient:

```powershell
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=phase9_report_package_audit
```

The aliases `workflow_mode=phase5_repeat_context`, `workflow_mode=phase5c_repeat_qc`, `workflow_mode=phase6_telomere`, `workflow_mode=phase7_rnaseq_readiness`, `workflow_mode=phase7b_quantification_strategy`, `workflow_mode=phase7c_expression_reference`, `workflow_mode=phase7d_rnaseq_quantification`, `workflow_mode=phase7e_expression_hardening`, `workflow_mode=phase8a_expression_integration`, `workflow_mode=phase8b_final_evidence_scoring`, and `workflow_mode=phase9_report_generation` remain available for stopping at earlier stages.

Use Snakemake only as a controlled single command, preferably outside native Windows or after confirming no previous Snakemake process is running. It requires the external OrthoFinder result tables to already be present or parsed, because OrthoFinder itself is run externally on this laptop. See `docs/runbook.md` for interruption recovery and WSL/Linux notes.

The June 2026 documentation/source cleanup rationale is recorded in `docs/maintainability_review.md`.

Run tests:

```bash
python -m pytest -q
```

## Expected Outputs

- `data/metadata/data_manifest.tsv`
- `data/metadata/species_manifest.tsv`
- `data/metadata/resource_status.tsv`
- `data/metadata/figshare_annotation_inventory.tsv`
- `results/qc/assembly_qc.tsv`
- `results/qc/protein_qc.tsv`
- `results/qc/reference_assembly_qc.tsv`
- `results/qc/reference_annotation_qc.tsv`
- `results/qc/reference_protein_qc.tsv`
- `results/claims/publication_claim_audit.tsv`
- `results/readiness/phase3_source_checks.tsv`
- `results/orthology/reference_gene_coordinates.tsv`
- `results/orthology/orthofinder_input_manifest.tsv`
- `results/orthology/orthofinder_preflight.tsv`
- `results/orthology/orthogroup_gene_counts_long.tsv` when `workflow_mode=orthology` completes
- `results/orthology/orthofinder_species_summary.tsv` when `workflow_mode=orthology` completes
- `results/orthology/candidate_copy_number.tsv`
- `results/orthology/candidate_duplication_audit.tsv`
- `results/evidence/phase3_integrated_evidence.tsv` when `workflow_mode=phase4d_consolidation` is used
- `results/evidence/phase4d_integrated_evidence.tsv` when `workflow_mode=phase5_repeat_context` is used
- `results/evidence/integrated_evidence.tsv`
- `results/evidence/phase4d_candidate_interpretation.tsv`
- `results/validation/candidate_isoform_audit.tsv`
- `results/validation/domain_check_targets.tsv`
- `results/validation/high_priority_rescue_targets.tsv`
- `results/validation/hmmer_pfam_preflight.tsv`
- `results/validation/candidate_domain_integrity.tsv`
- `results/domains/representative_candidate_proteins.pfam.domtblout`
- `data/interim/candidate_validation/representative_candidate_proteins.faa`
- `results/rescue/phase4_rescue_query_inventory.tsv`
- `results/rescue/phase4_forward_hits.tsv`
- `results/rescue/phase4_reciprocal_hits.tsv`
- `results/rescue/phase4_annotation_rescue_summary.tsv`
- `results/rescue/phase4_genome_rescue_plan.tsv`
- `data/interim/annotation_rescue/phase4_rescue_candidate_proteins.faa`
- `results/rescue/phase4b_miniprot_preflight.tsv`
- `data/interim/annotation_rescue/phase4b_miniprot_queries.faa`
- `data/interim/annotation_rescue/phase4b_miniprot_queries.tsv`
- `data/interim/annotation_rescue/phase4b_target_scaffolds.fna`
- `results/rescue/phase4b_target_regions.tsv`
- `results/rescue/phase4b_miniprot_raw.gff`
- `results/rescue/phase4b_miniprot.stderr.log`
- `results/rescue/phase4b_genome_alignment_hits.tsv`
- `results/rescue/phase4b_genome_validation_summary.tsv`
- `results/rescue/phase4c_locus_review.tsv`
- `results/rescue/phase4c_gene_review_summary.tsv`
- `results/rescue/tp53_targeted_forward_search_summary.tsv`
- `results/rescue/phase4e_locus_manual_review.tsv`
- `results/rescue/phase4e_gene_hardened_summary.tsv`
- `results/evidence/phase4e_hardened_evidence.tsv`
- `data/metadata/phase5_repeat_resource_status.tsv`
- `results/repeats/phase5_repeat_features.tsv`
- `results/repeats/phase5_candidate_locus_repeat_context.tsv`
- `results/repeats/phase5_gene_repeat_context_summary.tsv`
- `results/repeats/phase5c_repeatmasker_integrity.tsv`
- `results/repeats/phase5c_locus_repeat_qc.tsv`
- `results/repeats/phase5c_gene_repeat_qc_summary.tsv`
- `results/evidence/phase5_repeat_context_evidence.tsv`
- `results/telomere/phase6_telomeric_motif_scan.tsv`
- `results/telomere/phase6_scaffold_end_enrichment.tsv`
- `results/telomere/phase6_telomere_gene_audit.tsv`
- `data/metadata/rnaseq_manifest.tsv`
- `results/rnaseq/phase7_rnaseq_readiness.tsv`
- `results/rnaseq/phase7_candidate_expression_plan.tsv`
- `results/rnaseq/phase7b_reference_quantification_strategy.tsv`
- `results/rnaseq/phase7b_candidate_quantification_map.tsv`
- `results/rnaseq/phase7b_quantification_run_plan.tsv`
- `data/interim/rnaseq/phase7c/smic_tokyo.annotation_transcripts.fna`
- `data/interim/rnaseq/phase7c/smic_tokyo.tx2gene.tsv`
- `results/rnaseq/phase7c_expression_reference_qc.tsv`
- `results/rnaseq/phase7c_candidate_reference_validation.tsv`
- `results/rnaseq/phase7d_raw_read_intake.tsv`
- `results/rnaseq/phase7d_fastq_qc.tsv`
- `results/rnaseq/phase7d_salmon_preflight.tsv`
- `results/rnaseq/phase7d_salmon_quant_summary.tsv`
- `results/rnaseq/phase7d_candidate_expression_matrix.tsv`
- `results/rnaseq/phase7e_run_qc_review.tsv`
- `results/rnaseq/phase7e_candidate_expression_hardened.tsv`
- `results/rnaseq/phase7e_parameter_review.tsv`
- `results/evidence/phase8a_expression_support_evidence.tsv`
- `results/evidence/phase8a_expression_integration_audit.tsv`
- `results/evidence/phase8a_integrated_evidence.tsv`
- `results/evidence/phase8b_final_integrated_evidence.tsv`
- `results/evidence/phase8b_tier_audit.tsv`
- `results/evidence/phase8b_mechanism_summary.tsv`
- `results/reporting/phase9_figure_manifest.tsv`
- `results/reporting/phase9_key_findings.tsv`
- `reports/figures/data/phase9_evidence_tier_summary.tsv`
- `reports/figures/data/phase9_reporting_class_summary.tsv`
- `reports/figures/data/phase9_mechanism_evidence_matrix.tsv`
- `reports/figures/data/phase9_artifact_context_summary.tsv`
- `results/repeats/phase5b_repeat_annotation_preflight.tsv`
- `results/repeats/phase5b_repeat_annotation_plan.tsv`
- `data/metadata/phase5b_repeat_annotation_provenance.tsv`
- `data/metadata/phase5b_repeat_output_inventory.tsv`
- `reports/generated/preliminary_report.md`
- `reports/generated/phase5_repeat_context_report.md`
- `reports/generated/phase5c_repeat_context_qc_report.md`
- `reports/generated/phase6_telomere_report.md`
- `reports/generated/phase7_rnaseq_readiness_report.md`
- `reports/generated/phase7b_quantification_strategy_report.md`
- `reports/generated/phase7c_expression_reference_report.md`
- `reports/generated/phase7d_rnaseq_quantification_report.md`
- `reports/generated/phase7e_expression_hardening_report.md`
- `reports/generated/phase8a_expression_evidence_report.md`
- `reports/generated/phase8b_final_evidence_report.md`
- `reports/final/greenland_shark_longevity_phase9_report.md`
- `reports/final/phase9_report_package_audit.md`
- `reports/figures/phase9_evidence_tier_summary.svg`
- `reports/figures/phase9_reporting_class_summary.svg`
- `reports/figures/phase9_mechanism_evidence_matrix.svg`
- `reports/figures/phase9_artifact_context_summary.svg`

## Repository Layout

```text
config/      Workflow configuration and candidate panels
data/demo/   Tiny artificial fixtures for runnable workflow plumbing
docs/        Study design, manifest notes, claims register, hypotheses
src/         Small reusable Python modules and CLIs
tests/       Lightweight pytest coverage
workflow/    Snakemake entry point
```

## License And Citation

The repository code and documentation are released under the MIT License; see `LICENSE`.

If you use this workflow or adapt its reporting structure, cite the repository using `CITATION.cff`. Public biological datasets referenced by the workflow retain their original database, publication, and provider terms.
