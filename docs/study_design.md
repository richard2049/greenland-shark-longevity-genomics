# Study Design

## Research Question

What public genomic and transcriptomic evidence supports candidate mechanisms of extreme longevity and long-term tissue maintenance in the Greenland shark, and which signals remain uncertain because of annotation quality, assembly differences, missing metadata, limited sample size, or lack of functional validation?

## Dataset Inclusion Criteria

- Publicly accessible genome assembly, annotation, protein set, transcriptomic dataset, or associated publication.
- Stable accession, DOI, resource URL, or repository link.
- Sufficient metadata to determine species, resource type, and intended use.
- Large raw data are registered before download and processed only when a workflow phase explicitly supports them.

## Comparator Logic

Comparator species should be chosen for clear biological and technical reasons:

- Chondrichthyan relatives with public genomes and protein sets.
- Species used by published Greenland shark studies where accessions and methods are clear.
- Outgroups suitable for orthology, domain, or gene-family interpretation.
- Comparator choices must separate aging/longevity hypotheses from deep-sea, cold-adapted, retinal, immune, or assembly-quality confounders.

## Mechanism Categories

The MVP candidate panels cover DNA repair and genome stability, p53 pathway, chromatin, telomere biology, ferroptosis and iron handling, antioxidant response, mitochondrial maintenance, immune and inflammatory regulation, proteostasis and autophagy, and retinal preservation.

## Evidence Tiers

| Tier | Meaning | Minimum interpretation |
|---|---|---|
| Tier 1 | Robust computational evidence | Orthology-supported, domain-compatible, and not readily explained by isoforms, fragments, or annotation artifacts. Duplication claims require separable loci or comparable genomic support. |
| Tier 2 | Plausible but incomplete evidence | Some supporting evidence exists, but critical validation is missing. |
| Tier 3 | Exploratory/speculative | Useful for hypothesis generation only. Symbol-only matches and demo-only outputs belong here at most. |
| Artifact/uncertain | Do not claim biologically | Evidence is conflicted, artifact-prone, or insufficient for a biological statement. |

BUSCO and assembly contiguity are resource-quality evidence only. They are not biological evidence for longevity mechanisms.

## Interpretation Guardrails

- Do not use "activated" unless transcriptomic design supports activity-state inference.
- Do not use "inactivated" unless direct loss-of-function evidence is present.
- Do not equate "not annotated" with "absent".
- Do not equate copy-number expansion with functional advantage.
- Do not equate association with causation.
- Do not claim human translational relevance outside a clearly marked speculative section.

## MVP Limitations

The MVP still uses artificial demo data for default plumbing, but selected public reference files are now registered, downloaded or manually inventoried, QCed, and staged for Phase 1/Phase 3 analyses. The repository has parsed a real OrthoFinder orthogroup run generated externally in Docker, run HMMER/Pfam domain validation on mapped representative proteins, added targeted Phase 4 rescue/validation layers, added limited Phase 5 repeat-context checks around Phase 4e loci, implemented Phase 6 telomere sequence-context readiness, registered Phase 7a RNA-seq metadata/readiness, added a Phase 7b candidate-expression quantification strategy, implemented Phase 7c transcript-reference construction/input validation, and added Phase 7d local raw-read intake/QC plus guarded Salmon quantification. It still does not download raw RNA-seq automatically, run differential expression, run default RepeatMasker/RepeatModeler, run BUSCO, perform selection tests, or run full-genome miniprot rescue. Published findings remain external claims until the required orthology, domain, isoform, coordinate, genome-alignment, repeat-context, expression, and cross-resource checks support them.

## Publication Claim Audits

Published findings are tracked in `config/publication_claims.yaml` and audited into `results/claims/publication_claim_audit.tsv`. This keeps resource-quality observations separate from biological findings. For example, the Yang et al. PNAS 2026 linker-histone and `FTH1B` copy-number findings are registered as claims to reproduce, but remain `NOT_ASSESSED` until the workflow has suitable protein/annotation resources and orthology, domain, isoform, and locus validation.

## Phase 3 Input Readiness

The first real Phase 3 task is input discovery, not biological inference. The workflow records searched PNAS, bioRxiv, NCBI, University of Tokyo, Figshare, and indexed repository routes in `results/readiness/phase3_source_checks.tsv`. It then writes `results/orthology/reference_gene_coordinates.tsv` and `results/orthology/orthofinder_input_manifest.tsv`.

The PNAS data availability section declares a Figshare genome-annotation source, including protein sequences predicted from transcriptome data. The manually downloaded Figshare files are inventoried with checksums before being promoted to workflow inputs.

`complete.proteins.faa` is treated as a candidate genome-annotation protein FASTA for protein QC and OrthoFinder staging. `complete.genomic.gff` is parsed for gene coordinates. The transcriptome-derived TransDecoder peptide FASTA and small CDS subsets are explicitly excluded from genome-wide orthology and copy-number analysis unless a later phase defines a specific use for them. This prevents the workflow from treating mixed-scope protein sets as gene-family evidence.

## Comparator Intake

The first real comparator set is deliberately chondrichthyan-first:

- `Callorhinchus milii` provides a holocephalan comparator.
- `Scyliorhinus canicula` provides an annotated shark comparator.
- `Amblyraja radiata` provides a batoid/skate comparator.
- `Rhincodon typus` provides a large-shark comparator.

All four are represented by RefSeq genome-wide protein FASTA files. This choice keeps the first OrthoFinder input set close to the focal lineage and reduces the risk that broad vertebrate divergence dominates early candidate-gene interpretation. Distant bony vertebrate outgroups can be added later if needed for rooting, broader gene-family context, or specific candidate-gene checks.

## Orthology Execution

OrthoFinder is the selected first-pass orthology engine because it uses whole-proteome sequence similarity to infer orthogroups across multiple species, which is the appropriate starting point for gene-family copy-number screening. Symbol matching is not sufficient for this repository because symbols are inconsistent across species and annotations. A custom reciprocal-best-hit workflow is also avoided at this stage because it would under-handle many-to-many gene families and duplications, which are central to the scientific question.

The workflow includes an `orthology` mode that runs OrthoFinder and parses `Orthogroups.GeneCount.tsv` into long-format copy counts and species summaries when OrthoFinder is available on `PATH`. On this native Windows laptop, the practical execution route is Docker: OrthoFinder is run externally in a `staphb/orthofinder:2.5.5` container, then the repository parses that result directory and runs `workflow_mode=orthology_postprocess`. OrthoFinder execution remains separate from `reference_only` because it is heavier, depends on external binaries, and should fail fast if the environment is incomplete.

When OrthoFinder has been run externally, the `orthology_postprocess` mode maps curated candidates onto the parsed OrthoFinder results. The mapping uses exact GFF gene symbols and compact aliases only; descriptive product names are not used as primary evidence because they can merge paralogs or transfer annotation noise. Candidate copy number is reported as unique Greenland shark gene loci, while raw OrthoFinder protein counts are retained separately to expose isoform inflation. These rows still cannot become Tier 1 evidence without domain integrity, isoform filtering, locus-level validation, and cross-resource support.

## Candidate-Specific Validation Preparation

The immediate post-orthology validation layer performs tasks that are technically defensible without additional heavy databases:

- It audits candidate isoforms and selects one representative protein per mapped Greenland shark locus using the longest-protein isoform rule.
- It writes representative candidate proteins to `data/interim/candidate_validation/representative_candidate_proteins.faa` for downstream domain analysis.
- It records domain-check targets in `results/validation/domain_check_targets.tsv` but does not infer domain integrity.
- It records high-priority rescue targets in `results/validation/high_priority_rescue_targets.tsv`.

This design separates completed local checks from deferred validation. Domain integrity requires InterProScan or HMMER/Pfam with recorded database versions, e-values, and coverage. Annotation rescue requires validated query/reference proteins and local genome FASTA before miniprot or equivalent protein-to-genome alignment is meaningful.

## Candidate Domain Validation

The first implemented domain-validation path uses HMMER against Pfam-A on `data/interim/candidate_validation/representative_candidate_proteins.faa`. On the current native Windows environment this is run through `pyhmmer`, because command-line HMMER is not available through the current `win-64` conda channels. This is narrower than InterProScan by design. The current goal is to check whether representative candidate proteins have qualifying conserved Pfam domains with explicit e-values and domain coverage, not to perform broad functional annotation.

The workflow mode `domain_validation` writes:

- `results/validation/hmmer_pfam_preflight.tsv`
- `results/domains/representative_candidate_proteins.pfam.domtblout`
- `results/domains/representative_candidate_proteins.pfam.txt`
- `results/validation/candidate_domain_integrity.tsv`

Domain classification is conservative:

- `DOMAIN_SUPPORTED`: at least one Pfam hit passes the configured independent e-value and full-domain HMM coverage thresholds.
- `PARTIAL_DOMAIN`: a significant Pfam hit is detected, but full-domain coverage is incomplete.
- `NO_EXPECTED_DOMAIN_DETECTED`: no qualifying Pfam hit is detected under the configured thresholds.
- `NOT_ASSESSED`: the representative sequence or domain run is not available.

These categories describe domain evidence only. They do not establish protein function, pathway activity, gene duplication validity, gene absence, loss of function, or aging relevance. InterProScan can be added later as a broader independent annotation layer once the HMMER/Pfam result is stable.

## Phase 4 Targeted Annotation Rescue

Phase 4 currently focuses only on high-priority genes that exact-symbol mapping left as annotation uncertainty: `H1F0`, `FTH1B`, `TP53`, and `RAD51`. The implemented layer uses reciprocal protein similarity with HMMER `phmmer` through `pyhmmer`.

This choice is intentionally conservative and staged. Reciprocal protein similarity is appropriate for prioritizing candidate focal protein models when product names or symbols are incomplete, but it is not a replacement for orthology, domain integrity, separable-locus inspection, or protein-to-genome alignment.

The workflow mode `annotation_rescue` writes:

- `results/rescue/phase4_rescue_query_inventory.tsv`
- `results/rescue/phase4_forward_hits.tsv`
- `results/rescue/phase4_reciprocal_hits.tsv`
- `results/rescue/phase4_annotation_rescue_summary.tsv`
- `results/rescue/phase4_genome_rescue_plan.tsv`
- `data/interim/annotation_rescue/phase4_rescue_candidate_proteins.faa`

The current protein-level rescue summary identifies candidate focal protein models for `H1F0`, `FTH1B`, and `RAD51`. `TP53` remains unresolved by reciprocal protein rescue under the configured filters. None of these statuses support absence, inactivation, pathway activity, validated duplication, or longevity-mechanism claims.

## Phase 4b Genome-Level Validation

The workflow mode `genome_validation` implements the next gate for Phase 4: splice-aware protein-to-genome alignment with miniprot GFF output. This method is used because Phase 4b needs genomic coordinates, exon/CDS structure, and possible frameshift/stop flags. BLAST-like similarity is not sufficient for this validation layer because local similarity alone cannot classify a candidate gene model or distinguish separable genomic loci.

The mode writes:

- `results/rescue/phase4b_miniprot_preflight.tsv`
- `data/interim/annotation_rescue/phase4b_miniprot_queries.faa`
- `data/interim/annotation_rescue/phase4b_miniprot_queries.tsv`
- `data/interim/annotation_rescue/phase4b_target_scaffolds.fna`
- `results/rescue/phase4b_target_regions.tsv`
- `results/rescue/phase4b_miniprot_raw.gff`
- `results/rescue/phase4b_miniprot.stderr.log`
- `results/rescue/phase4b_genome_alignment_hits.tsv`
- `results/rescue/phase4b_genome_validation_summary.tsv`

On the current Windows laptop, miniprot is run through `tools/miniprot-docker.cmd` and the configured genome FASTA is the local Tokyo assembly file. The practical profile is `candidate_scaffolds`: the workflow extracts scaffolds implicated by Phase 4 protein-level rescue candidates and aligns the Phase 4b query proteins there. This avoids the memory cost of a full-genome Docker run, but it also means genes without candidate rescue coordinates, such as unresolved `TP53`, remain `not_assessable`.

The current Phase 4b summary reports separable high-coverage candidate-scaffold alignments for `H1F0`, `FTH1B`, and `RAD51`. These are candidate genome-alignment observations only. They do not validate duplication, function, pathway activity, adaptation, absence, inactivation, or aging relevance without manual exon/locus review, domain support, cross-resource checks, and conservative integration in Phase 8.

## Phase 4c/4d Locus Review And Consolidation

Phase 4c is a structured review layer over existing Phase 4b and Pfam outputs. It inspects overlap among candidate miniprot loci, local coordinates, Pfam support for rescue queries, and parsed miniprot disruption tags. `TP53` is handled separately with a bounded targeted search against scaffolds selected from Phase 4 forward-hit coordinates, because a broad chunked full-genome pass was not stable on the current laptop.

Phase 4d consolidates those validation tables into:

- `results/evidence/phase4d_candidate_interpretation.tsv`
- `results/evidence/integrated_evidence.tsv`
- `results/evidence/phase3_integrated_evidence.tsv`

The Phase 4d rule is intentionally deterministic: it joins already generated validation tables and applies conservative tiering rules. This is preferable to a new discovery algorithm at this stage because the task is interpretation and auditability, not additional candidate search. Current Phase 4d output treats `H1F0`, `FTH1B`, and `RAD51` as Tier 2 candidate-locus findings requiring manual review, and treats `TP53` as `Artifact/uncertain`. No Phase 4d row validates copy-number expansion, function, pathway activity, adaptation, or translational relevance.

## Phase 4e Manual Locus Hardening

Phase 4e turns the Phase 4c/4d candidate-locus review into a more explicit manual-review table. It does not run new discovery searches. Instead, it joins existing miniprot, Pfam, and targeted `TP53` outputs to the local Figshare annotation GFF and records:

- exact/focal annotation overlap,
- product-consistent gene-family overlap,
- candidate-locus spacing on the same scaffold,
- miniprot disruption flags,
- domain support,
- artifact flags,
- and conservative gene-level hardening status.

The current Phase 4e output writes:

- `results/rescue/phase4e_locus_manual_review.tsv`
- `results/rescue/phase4e_gene_hardened_summary.tsv`
- `results/evidence/phase4e_hardened_evidence.tsv`

Current interpretation remains conservative. `H1F0` and `RAD51` have candidate loci with miniprot, Pfam, and focal annotation support, but duplication status is not validated. `FTH1B` has a ferritin-family cluster with 12 domain-supported loci, 11 product-consistent annotation overlaps, 8 exact/focal annotation overlaps, and one miniprot disruption flag; this is Tier 2 candidate-locus evidence but high artifact risk. `TP53` has a product-consistent p53-family targeted alignment with a miniprot disruption flag and remains `Artifact/uncertain`.

## Phase 5 Limited Repeat Context

Phase 5 is currently limited to candidate-locus context from existing or externally generated repeat annotations. It intentionally does not run de novo repeat discovery inside the default workflow. This choice is methodologically conservative because repeat discovery is a heavier, parameter-sensitive workflow and should not be mixed into the high-priority candidate audit until the resource and interpretation layers are stable.

The implemented Phase 5 step checks the local Tokyo/PNAS Figshare inventory and configured repeat-annotation paths, normalizes RepeatMasker GenBank accessions to assembly-report scaffold names, filters repeat records to Phase 4e candidate windows, and intersects those intervals with the Phase 4e loci for `FTH1B`, `H1F0`, `RAD51`, and `TP53`.

Current Phase 5 output writes:

- `data/metadata/phase5_repeat_resource_status.tsv`
- `results/repeats/phase5_repeat_features.tsv`
- `results/repeats/phase5_candidate_locus_repeat_context.tsv`
- `results/repeats/phase5_gene_repeat_context_summary.tsv`
- `results/evidence/phase5_repeat_context_evidence.tsv`
- `reports/generated/phase5_repeat_context_report.md`

The current Phase 5 output imports the externally generated RepeatMasker GFF from Phase 5b and records local repeat context around the hardened candidate loci. `FTH1B`, `RAD51`, and `TP53` have direct repeat overlap in the bounded candidate-locus table, while `H1F0` has repeat annotations in the local window but no direct overlap. These observations are artifact/context evidence only and must not be interpreted as transposable-element mechanism, validated duplication, adaptation, pathway activity, or longevity evidence.

## Phase 5c Repeat-Context QC Hardening

Phase 5c hardens the Phase 5 repeat-context result without rerunning repeat discovery or masking. It streams the RepeatMasker `.out` file through the same Phase 4e candidate windows, compares direct-overlap calls with the GFF-derived Phase 5 table, records whether the RepeatMasker `.gff`, `.out`, and `.tbl` files are available, and keeps the `ProcessRepeats` signal-9 warning visible.

Current Phase 5c output writes:

- `results/repeats/phase5c_repeatmasker_integrity.tsv`
- `results/repeats/phase5c_locus_repeat_qc.tsv`
- `results/repeats/phase5c_gene_repeat_qc_summary.tsv`
- `reports/generated/phase5c_repeat_context_qc_report.md`

The `.out` parser supports the direct-overlap calls for `FTH1B`, `RAD51`, and `TP53`, while `H1F0` remains local-window repeat context without direct locus overlap. Because the `.tbl` summary is missing and the log records a signal-9 warning, the repeat outputs remain suitable for bounded candidate-locus artifact/context review only, not genome-wide repeat percentages or repeat-expansion claims.

## Phase 5b Optional De Novo Repeat Annotation

Phase 5b is the external route used to replace the earlier repeat-context `NOT_ASSESSED` status with a reproducible candidate-locus repeat-annotation table. It is optional and external by design. The repository writes preflight, command-plan, provenance, and output-inventory tables for a Docker/WSL/Linux RepeatModeler2 plus RepeatMasker workflow, but it does not run those heavy tools automatically.

The current Phase 5b planning outputs are:

- `results/repeats/phase5b_repeat_annotation_preflight.tsv`
- `results/repeats/phase5b_repeat_annotation_plan.tsv`
- `data/metadata/phase5b_repeat_annotation_provenance.tsv`
- `data/metadata/phase5b_repeat_output_inventory.tsv`

The preflight confirms that the local Tokyo genome FASTA is available and readable. The current external run produced a RepeatModeler family library and parseable RepeatMasker `.gff`/`.out` files, but no `.tbl` genome-wide summary. The RepeatMasker log includes a `ProcessRepeats` signal-9 warning, so the outputs are used for bounded candidate-locus repeat context rather than genome-wide repeat percentages or repeat-class claims.

This staged design keeps method-heavy repeat discovery separate from candidate-locus interpretation and prevents repeat overlap from being promoted to biological evidence without quality checks.

## Phase 6 Telomere Motif And Gene Readiness

Phase 6 implements a limited telomere-related analysis. It uses exact canonical motif counting for `TTAGGG` and `CCCTAA` in scaffold-end windows and adjacent internal control windows, plus a telomere/shelterin gene-readiness audit from the current candidate and integrated evidence tables.

The exact-motif approach is deliberately narrow. It is appropriate for detecting whether canonical motifs are enriched near scaffold ends, but it is not a telomere-length estimator and does not measure telomerase activity. More specialized telomere tools or long-read validation would be required before any telomere biology interpretation.

Current Phase 6 output writes:

- `results/telomere/phase6_telomeric_motif_scan.tsv`
- `results/telomere/phase6_scaffold_end_enrichment.tsv`
- `results/telomere/phase6_telomere_gene_audit.tsv`
- `reports/generated/phase6_telomere_report.md`

The current Tokyo assembly scan recorded 2,136 sequences, a 10 kb terminal window, and terminal motif enrichment relative to adjacent internal control windows (`end_enrichment_ratio = 1.61675`). The telomere-gene audit records current single-copy exact-symbol candidate mappings for `TERT`, `TERF1`, `TERF2`, `POT1`, `TINF2`, `ACD`, and `TERF2IP`. These are sequence-context and gene-readiness observations only. They must not be interpreted as telomere length, telomerase activity, rejuvenation, pathway activity, causation, or longevity-mechanism evidence.

## Phase 7a RNA-seq Metadata And Expression Readiness

Phase 7a registers run-level metadata for the public `PRJNA1246101` resource before any raw-read processing. The workflow uses NCBI SRA RunInfo metadata because this phase asks whether the public transcriptomic design is suitable for future candidate-expression support, not whether any gene is currently expressed.

Current Phase 7a output writes:

- `data/metadata/rnaseq_manifest.tsv`
- `results/rnaseq/phase7_rnaseq_readiness.tsv`
- `results/rnaseq/phase7_candidate_expression_plan.tsv`
- `reports/generated/phase7_rnaseq_readiness_report.md`

The current readiness table records three PolyA paired-end retinal RNA-seq runs (`SRR32965275`, `SRR32965277`, and `SRR32965276`) as metadata-ready for future candidate-expression audit, and excludes the paired-end WGS run (`SRR32965274`) from expression planning. All candidate-panel genes are staged for a future retina-only expression audit with `NOT_QUANTIFIED_PHASE7A_METADATA_ONLY` status. This is not expression support and must not be used for differential expression, pathway activity, whole-organism aging interpretation, or longevity-mechanism claims.

## Phase 7b Candidate Expression Quantification Strategy

Phase 7b converts Phase 7a readiness into a conservative future quantification plan. It maps candidate-panel genes to current annotation, CDS, isoform, copy-number, and Phase 4e hardening outputs, then classifies whether each candidate is ready for first-pass quantification or requires targeted locus/paralog review before expression-support interpretation.

The preferred future first-pass method is Salmon-style transcript quantification from an annotation-derived transcript FASTA with genome decoys and a transcript-to-gene map. This is chosen because the immediate public expression design is small, paired-end, and retina-specific, and because Salmon preserves transcript-to-gene ambiguity information useful for candidate-level summaries. Genome-aligned counting with STAR or HISAT2 plus featureCounts is reserved for ambiguous candidates where locus specificity and multi-mapping need direct inspection.

Phase 7b writes:

- `results/rnaseq/phase7b_reference_quantification_strategy.tsv`
- `results/rnaseq/phase7b_candidate_quantification_map.tsv`
- `results/rnaseq/phase7b_quantification_run_plan.tsv`
- `reports/generated/phase7b_quantification_strategy_report.md`

These are strategy/readiness outputs only. They do not show that any candidate gene is detected or expressed, and they do not support differential expression, pathway activity, tissue-general conclusions, causation, or longevity-mechanism claims.

## Phase 7c Expression Reference Construction

Phase 7c constructs an annotation-derived transcript FASTA and `tx2gene` table from the Tokyo genome FASTA and Figshare GFF. It then validates candidate gene IDs from Phase 7b against that reference before any raw-read quantification is attempted.

The implemented method extracts GFF3 `exon` intervals and streams the genome FASTA one scaffold at a time. This avoids the memory cost of loading the full assembly while keeping the construction deterministic and testable. CDS intervals are used only as a recorded fallback for transcript models without exon intervals.

Phase 7c writes:

- `data/interim/rnaseq/phase7c/smic_tokyo.annotation_transcripts.fna`
- `data/interim/rnaseq/phase7c/smic_tokyo.tx2gene.tsv`
- `results/rnaseq/phase7c_expression_reference_qc.tsv`
- `results/rnaseq/phase7c_candidate_reference_validation.tsv`
- `reports/generated/phase7c_expression_reference_report.md`

These outputs validate reference construction and candidate ID mapping only. They do not show expression, differential expression, pathway activity, or longevity relevance.

## Phase 7d Retina RNA-seq Intake/QC And Quantification

Phase 7d is a guarded execution layer for local retinal RNA-seq FASTQs. It does not fetch large SRA files automatically. Instead, it inspects whether expected paired FASTQs are present under the configured ignored raw-data directory, records file status and sizes, performs lightweight FASTQ sanity QC, checks Salmon executable/index readiness, and runs Salmon only when all required inputs are present.

The first-pass quantification method is Salmon against the Phase 7c annotation-derived transcript FASTA. This is appropriate for a small retina-only candidate-expression audit because Salmon is fast, uses probabilistic transcript assignment, and preserves transcript-level ambiguity for later `tx2gene` summarization. This is not a replacement for locus-aware genome alignment when paralogy, duplicated loci, repeats, or ambiguous gene families are central to interpretation. High-risk candidates such as `FTH1B`, `H1F0`, `RAD51`, and `TP53` still require targeted multi-mapping/locus review before expression support is used in Phase 8.

Phase 7d writes:

- `results/rnaseq/phase7d_raw_read_intake.tsv`
- `results/rnaseq/phase7d_fastq_qc.tsv`
- `results/rnaseq/phase7d_salmon_preflight.tsv`
- `results/rnaseq/phase7d_salmon_quant_summary.tsv`
- `results/rnaseq/phase7d_candidate_expression_matrix.tsv`
- `reports/generated/phase7d_rnaseq_quantification_report.md`

These outputs may support cautious retina-specific `detected` or `not detected under these conditions` language only after successful quantification and ambiguity review. They do not support differential expression, pathway activity, whole-organism aging interpretation, causation, functional advantage, or longevity-mechanism claims.

## Phase 7e/8a Expression Hardening And Evidence Integration

Phase 7e reviews completed Phase 7d quantification before expression support can be carried forward. It joins run-level mapping rates, lightweight FASTQ QC, Salmon index caveats, Phase 7c reference ambiguity, Phase 3/4 candidate validation, and Phase 5 repeat-context artifact risk. The current rules require consistent retina detection in at least two of three quantified runs for low-ambiguity candidates, while high-ambiguity candidates are deferred from positive expression support.

Phase 8a integrates the Phase 7e table into evidence scoring with a strict no-upgrade policy. It may append retina-specific expression support or caveats to a candidate row, but it preserves the existing biological evidence tier unless a later scoring step has explicit orthology, domain, locus, cross-resource, and expression-reference validation. This deterministic rule-based integration is preferred over a statistical score because the current data are retina-only and do not support a differential-expression design.

Phase 8a writes:

- `results/evidence/phase8a_expression_support_evidence.tsv`
- `results/evidence/phase8a_expression_integration_audit.tsv`
- `results/evidence/phase8a_integrated_evidence.tsv`
- `reports/generated/phase8a_expression_evidence_report.md`

## Phase 8b Final Evidence Scoring And Tier Audit

Phase 8b is the final pre-reporting evidence audit. It consolidates the Phase 8a evidence table with Phase 4e locus hardening, Phase 5c repeat-context QC, Phase 6 telomere-gene readiness, and Phase 7e/8a expression support. The goal is to classify each candidate as robust, plausible, exploratory, artifact-prone, or uncertain for Phase 9 reporting.

The method is deterministic rule-based auditing, not a numeric weighted score. This is preferred because the evidence sources are heterogeneous and incomplete: orthology, protein domains, miniprot/locus context, repeats, telomere motif context, and retina expression do not share a common statistical scale. Phase 8b can retain or conservatively downgrade evidence tiers, but it does not upgrade tiers. No current Phase 8b row reaches Tier 1.

Phase 8b writes:

- `results/evidence/phase8b_final_integrated_evidence.tsv`
- `results/evidence/phase8b_tier_audit.tsv`
- `results/evidence/phase8b_mechanism_summary.tsv`
- `reports/generated/phase8b_final_evidence_report.md`

## Tables And Figures

Machine-readable tables belong in the phase that creates the underlying evidence: Phase 1 for resource QC, Phase 3 for orthogroups/copy number, Phase 4 for validation/rescue, Phase 6 for telomere sequence context, Phase 7 for expression metadata and later expression support, and Phase 8 for integrated evidence scoring. Publication-quality figures should be concentrated in Phase 9 after the underlying tables are stable. Earlier phases can include diagnostic plots, but the report-level figures should be generated from versioned result tables so that interpretation remains traceable.

## Phase 9 Report And Figure Generation

Phase 9 is a deterministic reporting layer over Phase 8b. It uses `results/evidence/phase8b_tier_audit.tsv` and `results/evidence/phase8b_mechanism_summary.tsv` as primary inputs, with `results/evidence/phase8b_final_integrated_evidence.tsv` retained for row-level traceability.

The repository generates categorical SVG figures and matching TSV figure-data tables instead of statistical plots or weighted scores. This is deliberate: Phase 8b contains evidence classes, artifact-risk classes, and reporting-readiness labels, not a replicate-level experimental design. Categorical bar charts and a mechanism evidence matrix make the current state interpretable while avoiding false precision.

Phase 9 figures follow conservative computational-biology reporting practices: use versioned source tables, keep plotted values machine-readable, avoid decorative or 3D encodings, use a colorblind-safe palette, label count axes directly, include source/provenance captions, and place biological interpretation guardrails on each figure. Context counters for expression, repeat overlap, and locus review are visually separated from evidence-tier bars so artifact/context information cannot be mistaken for stronger mechanism support.

Phase 9 writes:

- `reports/final/greenland_shark_longevity_phase9_report.md`
- `results/reporting/phase9_figure_manifest.tsv`
- `results/reporting/phase9_key_findings.tsv`
- `reports/figures/data/phase9_evidence_tier_summary.tsv`
- `reports/figures/data/phase9_reporting_class_summary.tsv`
- `reports/figures/data/phase9_mechanism_evidence_matrix.tsv`
- `reports/figures/data/phase9_artifact_context_summary.tsv`
- `reports/figures/phase9_evidence_tier_summary.svg`
- `reports/figures/phase9_reporting_class_summary.svg`
- `reports/figures/phase9_mechanism_evidence_matrix.svg`
- `reports/figures/phase9_artifact_context_summary.svg`

Phase 9 must not change evidence tiers, infer pathway state, validate duplication, infer telomere length, generalize retina expression to whole organism, infer causation, or make human translational claims. Its role is to answer what is robust under current criteria, what is plausible but incomplete, what is artifact-prone, what is uncertain, and what requires follow-up.

After Phase 9, the repository can run a report-package audit. This is a release-quality check over traceability, figure metadata/provenance, language guardrails, `.gitignore` behavior, and public-repository packaging risks. It is not a new analysis phase and does not change biological evidence tiers. The current audit records that TSV-backed SVGs are suitable for this categorical evidence-audit stage, while optional matplotlib/R exports may be useful later for manuscript-specific formatting.
