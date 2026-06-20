# Phase 4 Annotation Rescue

Phase 4 targets only high-priority candidate genes that remain unresolved because exact-symbol mapping did not identify a Greenland shark locus. This is an annotation-rescue step, not a biological claim.

## Scope

Current targets are read from `results/validation/high_priority_rescue_targets.tsv` and restricted by `config/config.yaml`:

- `H1F0`
- `FTH1B`
- `TP53`
- `RAD51`

All four entered Phase 4 with `ANNOTATION_UNCERTAINTY_NO_EXACT_SYMBOL_MATCH`. That status must not be interpreted as gene absence.

## Method

The first implemented rescue layer uses `pyhmmer`/HMMER `phmmer` reciprocal protein similarity:

1. Select comparator query proteins from staged OrthoFinder input proteomes using conservative gene-specific product-description patterns.
2. Search selected comparator query proteins against the focal Greenland shark protein FASTA.
3. Retain forward hits that pass e-value and coverage filters.
4. Search each passing focal hit back against the source comparator proteome.
5. Classify reciprocal top hits as protein-level rescue candidates only when the reciprocal hit returns to a selected query or gene-pattern-compatible query.

This protein-level approach is used before genome alignment because local staged protein sets provide a fast, inspectable filter for unresolved candidates. Protein-level rescue is sufficient to prioritize candidates for miniprot; it is not sufficient to validate a gene model or duplication.

## Outputs

- `results/rescue/phase4_rescue_query_inventory.tsv`: comparator query proteins selected or excluded by query cap.
- `results/rescue/phase4_forward_hits.tsv`: forward `phmmer` hits into the focal protein FASTA.
- `results/rescue/phase4_reciprocal_hits.tsv`: reciprocal top-hit status for passing forward hits.
- `results/rescue/phase4_annotation_rescue_summary.tsv`: per-gene conservative rescue status.
- `results/rescue/phase4_genome_rescue_plan.tsv`: genome-alignment readiness and blockers.
- `data/interim/annotation_rescue/phase4_rescue_candidate_proteins.faa`: focal candidate protein sequences for later domain and locus review.

## Current Result

The current run produced protein-level rescue candidates for `H1F0`, `FTH1B`, and `RAD51`. `TP53` had forward protein-similarity hits but no reciprocal protein-level rescue candidate under the configured filters.

These results should be interpreted as follows:

- `H1F0`: candidate focal proteins require domain/region review and locus inspection.
- `FTH1B`: candidate focal proteins require especially careful ferritin-family review because ferritin annotations are duplication- and isoform-sensitive.
- `RAD51`: candidate focal proteins include isoforms at candidate loci and require isoform/locus review.
- `TP53`: unresolved after reciprocal protein rescue; this does not support absence, inactivation, or loss.

## Guardrails

- Do not report `activated`, `inactivated`, `absent`, causal, adaptive, or validated-duplication language from these tables.
- Do not treat product names alone as orthology.
- Do not treat reciprocal protein rescue as genome-level annotation validation.
- Do not infer loss from failure to rescue.

## Phase 4b Genome-Validation Gate

The repository now implements a gated genome-level validation layer in `workflow_mode=genome_validation`.

The selected method is miniprot GFF output. This is preferred over BLAST-like local similarity for Phase 4b because the scientific question is whether candidate proteins can be aligned to interpretable genomic loci with exon structure and possible frameshift/stop signals. BLAST can identify local similarity, but it does not provide the same splice-aware protein-to-genome model needed for candidate gene-model validation.

Phase 4b writes:

- `results/rescue/phase4b_miniprot_preflight.tsv`: required input and executable checks.
- `data/interim/annotation_rescue/phase4b_miniprot_queries.faa`: focal rescue candidates plus comparator queries for unresolved genes.
- `data/interim/annotation_rescue/phase4b_miniprot_queries.tsv`: query provenance table.
- `data/interim/annotation_rescue/phase4b_target_scaffolds.fna`: candidate-scaffold FASTA used by the current laptop profile.
- `results/rescue/phase4b_target_regions.tsv`: scaffold extraction provenance.
- `results/rescue/phase4b_miniprot_raw.gff`: raw miniprot GFF, or a documented not-run GFF stub.
- `results/rescue/phase4b_miniprot.stderr.log`: miniprot stderr, or a not-run explanation.
- `results/rescue/phase4b_genome_alignment_hits.tsv`: parsed alignment hits when miniprot runs.
- `results/rescue/phase4b_genome_validation_summary.tsv`: conservative per-gene genome-validation status.

Current Phase 4b status is based on a limited `candidate_scaffolds` run against the local Tokyo genome FASTA with `miniprot` provided by `tools/miniprot-docker.cmd`. The workflow prepares 15 query proteins: focal rescue-candidate proteins for `H1F0`, `FTH1B`, and `RAD51`, plus the selected comparator query for unresolved `TP53`.

The current summary detects separable high-coverage candidate-scaffold alignments for `H1F0`, `FTH1B`, and `RAD51`. `TP53` remains `not_assessable` in this run because it has no protein-level rescued candidate region and the workflow did not run a full-genome search for it.

Phase 4b classifies results conservatively as `intact_candidate`, `duplicated_candidate`, `possible_fragment`, `possible_pseudogene`, `annotation_uncertainty`, or `not_assessable`. These labels remain candidate validation statuses, not final biological claims. Miniprot terminal `stop_codon` features are treated as normal annotations; potential disruption is parsed from miniprot PAF tags such as `fs:i` and `st:i` when present.
