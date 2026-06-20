# Data Manifest Notes

The workflow writes machine-readable manifest files under `data/metadata/`. This document records the initial manually reviewed public resources and how they should be used.

| Resource ID | Species | Accession | Source | File type | Status | Size if known | Usage notes | Publication/resource |
|---|---|---|---|---|---|---:|---|---|
| `SMIC_FLI_GENOME_2025` | *Somniosus microcephalus* | `GCA_052327205.1_FLI_Smic_1.0`; `PRJNA1236315`; `SAMN47384478` | NCBI/UCSC GenArk/FLI | Genome assembly/browser resource | `REGISTERED` | `TODO` | Manifest only in MVP. Use for future assembly/protein/annotation QC if downloadable files are selected. | Sahm et al. bioRxiv preprint and FLI/UCSC genome browser |
| `SMIC_TOKYO_GENOME_2025` | *Somniosus microcephalus* | `PRJNA1218902`; `GCA_056099535.1`; `JBLTJD000000000`; `SAMN46535397` | NCBI BioProject | Genome assembly and WGS/SRA resource | `REGISTERED` | `383 Gbases`; `0.16 Tbytes` BioProject data volume | Manifest only in MVP. Useful for cross-resource reproducibility checks after local download policy is defined. | University of Tokyo BioProject |
| `SMIC_RETINA_PRJNA1246101_2026` | *Somniosus microcephalus* and comparator sharks | `PRJNA1246101`; coding sequences `PV442159-PV442194` | NCBI BioProject/SRA; Nature Communications | Raw genomic and retinal RNA-seq reads; coding sequences | `REGISTERED` | `264 Gbases`; `85456 Mbytes` BioProject data volume | Run-level RNA-seq metadata are registered in Phase 7a. Phase 7d can process locally supplied retinal FASTQs, but raw reads are not downloaded automatically. | Fogg et al. Nature Communications 2026 |

## Phase 1 Reference Intake

The file-level intake table is generated at `data/metadata/reference_file_inventory.tsv`.

As of the 2026-05-27 intake pass:

- The exact NCBI genome FASTA and GenBank flatfile URLs are registered for `GCA_052327205.1_FLI_Smic_1.0` and `GCA_056099535.1_ASM5609953v1`, but these gigabyte-scale files are not downloaded by default.
- Small NCBI assembly statistics, assembly report, feature-count, and checksum files are downloaded under `data/raw/references/`.
- The registered NCBI assembly directories do not list separate protein FASTA or genomic GFF files. This is an annotation/resource-availability limitation, not evidence for gene loss.
- The `PV442159-PV442194` coding-sequence accessions from `PRJNA1246101` are downloaded as nucleotide FASTA and CDS-derived amino-acid FASTA for small-file QC. These files are not a genome-wide protein set.
- Phase 1 QC from `assembly_stats.txt` is written to `results/qc/reference_assembly_qc.tsv`.
- Protein sequence QC is `NOT_ASSESSED` for the two assembly packages, and is run only for the `PRJNA1246101` CDS-derived amino-acid FASTA in `results/qc/reference_protein_qc.tsv`.

Publication claim audits are generated at `results/claims/publication_claim_audit.tsv`. The PNAS/Yang et al. 2026 claims are registered there as claims requiring reproduction, not as accepted repository findings.

## Phase 3 Input Discovery

The file `results/readiness/phase3_source_checks.tsv` records the 2026-05-27 search for genome-wide protein FASTA, GFF/GTF annotation, gene models, and orthogroup/copy-number supplementary tables. Searched routes include:

- PNAS final article DOI `10.1073/pnas.2601272123` and candidate PNAS supplement endpoints.
- bioRxiv preprint DOI `10.1101/2025.02.19.638963` through the bioRxiv API and registered JATS XML endpoint.
- NCBI Datasets reports and FTP package inventories for `GCA_056099535.1` and `GCA_052327205.1`.
- University of Tokyo author-institution release.
- Indexed GitHub/web searches for associated repositories.

Manual inspection of the PNAS data availability section identified a Figshare source for genome annotation, including protein sequences predicted from transcriptome data: `https://figshare.com/s/4f1adabbc84fbf5a72e0`. This source is registered in `config/config.yaml` and `results/readiness/phase3_source_checks.tsv`.

The manually downloaded Figshare archive and extracted files are inventoried in `data/metadata/figshare_annotation_inventory.tsv` with file sizes, MD5 checksums, SHA-256 checksums, and FASTA sequence counts where applicable. The current files are:

- `complete.genomic.gff`: genome annotation GFF used for coordinate parsing.
- `complete.proteins.faa`: candidate genome-annotation protein FASTA used for real protein QC and OrthoFinder input staging.
- `complete.cds.fna`: coding-sequence FASTA, not used directly for OrthoFinder.
- `trinity_out_dir.Trinity.fasta.transdecoder.pep`: transcriptome-derived peptide FASTA, inventoried but kept separate from genome-wide copy-number analysis.
- `greenland shark annotation.zip`: original downloaded archive retained for provenance.

Coordinate and protein QC outputs are resource-quality/input-readiness observations only. They must not be interpreted as biological validation of any candidate mechanism without orthology, domain, isoform, and locus-level checks.

## Phase 3 Comparator Proteomes

The first comparator intake downloads selected RefSeq protein FASTA files for:

| Resource ID | Species | Accession | Role |
|---|---|---|---|
| `CMIL_REFSEQ_2021` | *Callorhinchus milii* | `GCF_018977255.1` | Holocephalan chondrichthyan comparator |
| `SCAN_REFSEQ_2026` | *Scyliorhinus canicula* | `GCF_902713615.2` | Shark comparator |
| `ARAD_REFSEQ_2026` | *Amblyraja radiata* | `GCF_010909765.2` | Batoid/skate comparator |
| `RTYP_REFSEQ_2022` | *Rhincodon typus* | `GCF_021869965.1` | Large-shark comparator |

These files are registered in `config/config.yaml`, downloaded into `data/raw/references/`, checksummed in `data/metadata/reference_file_inventory.tsv`, QCed in `results/qc/reference_protein_qc.tsv`, and decompressed/staged as plain FASTA files in `data/interim/orthofinder_input/`. This completes input readiness for a small chondrichthyan OrthoFinder run.

## Phase 3 OrthoFinder Run Provenance

The current orthogroups were generated outside Snakemake with Docker because native Windows conda could not provide OrthoFinder. The run recorded in `results/orthofinder_og_20260528_153948/Results_May28/Log.txt` used:

- Container image tag: `staphb/orthofinder:2.5.5`
- OrthoFinder version: `2.5.5`
- Command: `orthofinder -f data/interim/orthofinder_input -o results/orthofinder_og_20260528_153948 -S diamond -og -t 8 -a 8`
- Start time in OrthoFinder log: `2026-05-28 13:39:52`
- Completion time in OrthoFinder log: `2026-05-28 14:53:41`
- Result directory: `results/orthofinder_og_20260528_153948/Results_May28/`

The `-og` option was used intentionally to stop after orthogroup inference. This is sufficient for the current candidate copy-number screening step and avoids unnecessary gene-tree/species-tree computation in the MVP.

`results/orthology/orthofinder_preflight.tsv` records whether the staged inputs and OrthoFinder executable are available. In the current Windows laptop workflow, OrthoFinder was run externally with Docker and then parsed into repository-standard outputs. When `workflow_mode=orthology` runs in WSL/Linux, or when an external Docker run is parsed, outputs are written to:

- `results/orthology/orthogroup_gene_counts_long.tsv`
- `results/orthology/orthofinder_species_summary.tsv`

The Docker/external-run continuation is:

```powershell
python -m greenland_shark_longevity.orthofinder_workflow parse `
  --results-dir results/orthofinder_og_YYYYMMDD_HHMMSS `
  --manifest results/orthology/orthofinder_input_manifest.tsv `
  --gene-count-output results/orthology/orthogroup_gene_counts_long.tsv `
  --species-summary-output results/orthology/orthofinder_species_summary.tsv

snakemake --snakefile workflow/Snakefile --cores 1 --config workflow_mode=orthology_postprocess
```

Those tables are still computational evidence inputs, not final biological claims. Candidate copy-number interpretation requires isoform filtering, domain checks, and coordinate validation. The current postprocessing reports unique Greenland shark loci for candidate copy number and retains raw OrthoFinder protein counts separately to expose possible isoform inflation.

## Candidate Validation Prepared Outputs

The candidate-specific validation preparation step writes:

- `results/validation/candidate_isoform_audit.tsv`: mapped candidate loci, protein isoforms, protein lengths, selected representative isoform, and coordinate context.
- `data/interim/candidate_validation/representative_candidate_proteins.faa`: representative candidate protein sequences for downstream domain tools.
- `results/validation/domain_check_targets.tsv`: proteins ready for InterProScan or HMMER/Pfam domain validation. Domain integrity remains `NOT_ASSESSED` until an explicit domain tool/database run is performed.
- `results/validation/high_priority_rescue_targets.tsv`: unresolved high-priority candidates requiring reciprocal similarity search and protein-to-genome rescue.

The current high-priority rescue set is `H1F0`, `FTH1B`, `TP53`, and `RAD51`. These are rescue targets because exact-symbol OrthoFinder mapping did not resolve them; they must not be interpreted as absent.

## Domain Validation Inputs

The HMMER/Pfam validation path expects a local Pfam-A HMM database at the configured path:

- `data/raw/references/PFAM/Pfam-A.hmm`

This file is not downloaded by the default workflow because it is a large external database with its own release cycle. The current local installation is Pfam release `38.1`, retrieved on `2026-05-29` from `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`. The compressed SHA-256 is `D3D30C8E6801BFEDECF783408ECC98916F8F1DDA8974C6E51036FCBDD765F591`; the decompressed HMM SHA-256 is `A78A42D6FAF265B6BFCA59E8F062D06FAE6083CE2C6E335D7B381F20B82B7903`.

The current native Windows workflow uses the `pyhmmer` backend, so `hmmpress` index files are not required. If the backend is changed to command-line `hmmscan`, record the generated `hmmpress` index files too. The workflow writes environment/resource readiness to `results/validation/hmmer_pfam_preflight.tsv`.

Domain-validation outputs are:

- `results/domains/representative_candidate_proteins.pfam.domtblout`
- `results/domains/representative_candidate_proteins.pfam.txt`
- `results/validation/candidate_domain_integrity.tsv`

These outputs are protein-domain evidence, not biological mechanism evidence by themselves.

## Phase 4 Annotation Rescue Outputs

The targeted Phase 4 rescue workflow uses already staged protein resources and does not download new public data. It writes:

- `results/rescue/phase4_rescue_query_inventory.tsv`
- `results/rescue/phase4_forward_hits.tsv`
- `results/rescue/phase4_reciprocal_hits.tsv`
- `results/rescue/phase4_annotation_rescue_summary.tsv`
- `results/rescue/phase4_genome_rescue_plan.tsv`
- `data/interim/annotation_rescue/phase4_rescue_candidate_proteins.faa`

The configured local genome FASTA for genome-level rescue is now the Tokyo assembly FASTA at `data/raw/references/SMIC_TOKYO_GENOME_2025/GCA_056099535.1_ASM5609953v1_genomic.fna.gz`. The current Windows execution uses `tools/miniprot-docker.cmd` as a containerized miniprot wrapper. Phase 4 protein-level rescue still remains a prioritization layer: it does not validate gene loss, pseudogenization, duplication, or function without Phase 4b locus review and later evidence integration.

## Phase 4b Genome Validation Manifest Notes

The `genome_validation` mode does not download new data. It prepares query files, records the local genome/miniprot preflight, extracts target scaffolds when `phase4b_genome_validation.genome_scope` is `candidate_scaffolds`, and runs miniprot when the required inputs are present:

- `data/interim/annotation_rescue/phase4b_miniprot_queries.faa`
- `data/interim/annotation_rescue/phase4b_miniprot_queries.tsv`
- `data/interim/annotation_rescue/phase4b_target_scaffolds.fna`
- `results/rescue/phase4b_target_regions.tsv`
- `results/rescue/phase4b_miniprot_preflight.tsv`
- `results/rescue/phase4b_miniprot_raw.gff`
- `results/rescue/phase4b_miniprot.stderr.log`
- `results/rescue/phase4b_genome_alignment_hits.tsv`
- `results/rescue/phase4b_genome_validation_summary.tsv`

The current preflight reports 15 query proteins ready for alignment and uses the candidate-scaffold FASTA as the actual miniprot target. The target scaffold extraction currently includes the scaffolds carrying Phase 4 rescue candidates for `H1F0`, `FTH1B`, and `RAD51`; `TP53` has no candidate region in this limited-scope run and remains `not_assessable`.

## Phase 4c/4d Review And Evidence Outputs

Phase 4c and Phase 4d do not download new data. They consolidate existing Phase 4 rescue, Pfam, and miniprot outputs:

- `results/rescue/phase4c_rescue_domain_targets.tsv`
- `results/rescue/phase4c_rescue_domain_integrity.tsv`
- `results/rescue/phase4c_locus_review.tsv`
- `results/rescue/phase4c_gene_review_summary.tsv`
- `results/rescue/tp53_forward_hit_target_regions.tsv`
- `results/rescue/tp53_targeted_forward_alignment_hits.tsv`
- `results/rescue/tp53_targeted_forward_search_summary.tsv`
- `results/evidence/phase3_integrated_evidence.tsv`
- `results/evidence/phase4d_candidate_interpretation.tsv`
- `results/evidence/integrated_evidence.tsv`

`phase3_integrated_evidence.tsv` is the pre-consolidation evidence snapshot from the candidate duplication audit. `phase4d_candidate_interpretation.tsv` records high-priority candidate interpretations from Phase 4 validation tables. `integrated_evidence.tsv` is updated with conservative Phase 4d rows for `H1F0`, `FTH1B`, `RAD51`, and `TP53`.

## Phase 4e Manual Locus Hardening Outputs

Phase 4e does not download new data. It uses:

- `results/rescue/phase4c_locus_review.tsv`
- `results/rescue/phase4c_rescue_domain_integrity.tsv`
- `results/rescue/tp53_targeted_forward_alignment_hits.tsv`
- `results/rescue/tp53_forward_hit_target_regions.tsv`
- `data/raw/references/SMIC_TOKYO_GENOME_2025/figshare_annotation/greenland shark annotation/complete.genomic.gff`

It writes:

- `results/rescue/phase4e_locus_manual_review.tsv`
- `results/rescue/phase4e_gene_hardened_summary.tsv`
- `results/evidence/phase4e_hardened_evidence.tsv`

These tables record annotation overlap and artifact-risk information for manual review. They do not validate duplication, gene function, gene state, pathway activity, adaptation, or translational relevance.

## Phase 5 Limited Repeat Context Outputs

Phase 5 does not download new data and does not run RepeatMasker or RepeatModeler. It checks only existing local files from the Tokyo/PNAS Figshare package and any explicit repeat-annotation paths configured in `config/config.yaml`.

The current local inspection found:

- `complete.genomic.gff` was scanned for repeat-like GFF feature types and repeat attributes; no parseable repeat features were detected in that source.
- `complete.proteins.faa`, `complete.cds.fna`, and the TransDecoder peptide FASTA are not repeat-annotation candidates.
- The original Figshare archive is retained for provenance, while extracted files are inspected separately.
- The externally generated Phase 5b RepeatMasker GFF is now registered in `config/config.yaml` as the repeat-annotation source used for bounded Phase 5 candidate-locus context.

The Phase 5 outputs are:

- `data/metadata/phase5_repeat_resource_status.tsv`: file-by-file repeat-annotation availability and inspection status.
- `results/repeats/phase5_repeat_features.tsv`: parsed repeat intervals from the Phase 5b RepeatMasker GFF, filtered to Phase 4e candidate windows.
- `results/repeats/phase5_candidate_locus_repeat_context.tsv`: repeat-context status for each Phase 4e locus.
- `results/repeats/phase5_gene_repeat_context_summary.tsv`: gene-level repeat-context summary for `FTH1B`, `H1F0`, `RAD51`, and `TP53`.
- `results/repeats/phase5c_repeatmasker_integrity.tsv`: RepeatMasker file presence, checksum provenance, missing `.tbl`, and log-warning status.
- `results/repeats/phase5c_locus_repeat_qc.tsv`: locus-level comparison between the GFF-derived repeat context and a window-filtered RepeatMasker `.out` parse.
- `results/repeats/phase5c_gene_repeat_qc_summary.tsv`: gene-level Phase 5c QC summary used by the integrated evidence layer.

Current Phase 5 statuses should be interpreted as resource availability and artifact-context readiness. They do not support transposable-element mechanism claims, repeat-mediated duplication claims, or aging/longevity interpretation.

## Phase 6 Telomere Readiness Manifest Notes

Phase 6 uses the local Tokyo genome FASTA and current candidate/evidence tables. It does not download new public data and does not use telomere-specific wet-lab or long-read validation methods.

The Phase 6 outputs are:

- `results/telomere/phase6_telomeric_motif_scan.tsv`: exact `TTAGGG`/`CCCTAA` counts in scaffold-end and adjacent internal control windows.
- `results/telomere/phase6_scaffold_end_enrichment.tsv`: assembly-level motif-density and scaffold-end enrichment summary.
- `results/telomere/phase6_telomere_gene_audit.tsv`: telomere/shelterin candidate readiness from current orthology and integrated evidence tables.
- `reports/generated/phase6_telomere_report.md`: human-readable Phase 6 summary.

Current Phase 6 statuses should be interpreted as sequence-context and gene-readiness observations only. They do not support telomere length, telomerase activity, rejuvenation, pathway activity, gene-state, or aging/longevity mechanism claims.

## Phase 7a RNA-seq Readiness Manifest Notes

Phase 7a registers `PRJNA1246101` run-level metadata from NCBI SRA RunInfo. It does not download raw reads or quantify expression.

The Phase 7a outputs are:

- `data/metadata/rnaseq_manifest.tsv`: exact run-level metadata, accessions, retrieval date, source URL, expected local raw-read paths, and download status.
- `results/rnaseq/phase7_rnaseq_readiness.tsv`: run-level expression-readiness classification.
- `results/rnaseq/phase7_candidate_expression_plan.tsv`: candidate-panel expression-audit plan with all genes marked `NOT_QUANTIFIED_PHASE7A_METADATA_ONLY`.
- `reports/generated/phase7_rnaseq_readiness_report.md`: human-readable metadata/readiness summary.

Current registered runs:

| Run | Experiment | BioSample | Library strategy | Tissue | Phase 7a use |
|---|---|---|---|---|---|
| `SRR32965275` | `SRX28238856` | `SAMN47762445` | RNA-Seq, PolyA, paired | retina | Future candidate expression audit |
| `SRR32965277` | `SRX28238854` | `SAMN47762443` | RNA-Seq, PolyA, paired | retina | Future candidate expression audit |
| `SRR32965276` | `SRX28238855` | `SAMN47762444` | RNA-Seq, PolyA, paired | retina | Future candidate expression audit |
| `SRR32965274` | `SRX28238857` | `SAMN47762446` | WGS, paired | retina metadata label | Retained for provenance; excluded from expression planning |

These tables are metadata/readiness evidence only. They do not support expression detection, differential expression, pathway activity, whole-organism aging interpretation, or longevity-mechanism claims.

## Phase 7b Candidate Quantification Strategy Manifest Notes

Phase 7b does not download raw reads or quantify expression. It uses existing local and generated metadata to prepare a candidate-to-reference mapping and future quantification plan:

- `data/metadata/rnaseq_manifest.tsv`
- `results/rnaseq/phase7_candidate_expression_plan.tsv`
- `results/orthology/candidate_copy_number.tsv`
- `results/validation/candidate_isoform_audit.tsv`
- `results/rescue/phase4e_gene_hardened_summary.tsv`
- `results/rescue/phase4e_locus_manual_review.tsv`
- `data/raw/references/SMIC_TOKYO_GENOME_2025/figshare_annotation/greenland shark annotation/complete.cds.fna`

The Phase 7b outputs are:

- `results/rnaseq/phase7b_reference_quantification_strategy.tsv`: preferred, fallback, and validation quantification routes with required inputs and limitations.
- `results/rnaseq/phase7b_candidate_quantification_map.tsv`: candidate-panel mapping to annotation/CDS IDs, isoform/paralog caveats, and readiness for future expression summarization.
- `results/rnaseq/phase7b_quantification_run_plan.tsv`: planned raw-read download, QC, transcript reference generation, Salmon indexing, quantification, and candidate summarization steps.
- `reports/generated/phase7b_quantification_strategy_report.md`: human-readable strategy summary.

These outputs are quantification-planning evidence only. They do not indicate expression detection, differential expression, pathway activity, whole-organism aging interpretation, or longevity-mechanism support.

## Phase 7c Expression Reference Manifest Notes

Phase 7c constructs a local expression reference from files already registered in the repository:

- `data/raw/references/SMIC_TOKYO_GENOME_2025/GCA_056099535.1_ASM5609953v1_genomic.fna.gz`
- `data/raw/references/SMIC_TOKYO_GENOME_2025/figshare_annotation/greenland shark annotation/complete.genomic.gff`
- `results/rnaseq/phase7b_candidate_quantification_map.tsv`

The Phase 7c outputs are:

- `data/interim/rnaseq/phase7c/smic_tokyo.annotation_transcripts.fna`: transcript sequences reconstructed from GFF exon intervals and the Tokyo genome FASTA.
- `data/interim/rnaseq/phase7c/smic_tokyo.tx2gene.tsv`: transcript-to-gene mapping used for future Salmon summarization.
- `results/rnaseq/phase7c_expression_reference_qc.tsv`: transcript-reference construction QC.
- `results/rnaseq/phase7c_candidate_reference_validation.tsv`: candidate-panel validation against the constructed reference.
- `reports/generated/phase7c_expression_reference_report.md`: human-readable reference-construction summary.

These outputs are intermediate reference files and validation tables, not expression evidence. Raw reads remain undownloaded unless a later phase explicitly runs read intake and QC.

## Phase 7d Raw Read Intake/QC And Salmon Quantification Manifest Notes

Phase 7d uses the Phase 7a manifest, Phase 7a readiness table, and Phase 7c expression reference. It does not download raw reads automatically. Expected paired FASTQs should be placed or linked under:

- `data/raw/rnaseq/SMIC_RETINA_PRJNA1246101_2026/`

The expected naming patterns are configured in `config/config.yaml` under `phase7d_rnaseq_quantification.expected_fastq_patterns`. Current output paths are:

- `results/rnaseq/phase7d_raw_read_intake.tsv`: per-run local FASTQ status, expected paths, file sizes, and required action.
- `results/rnaseq/phase7d_fastq_qc.tsv`: lightweight FASTQ record-count, read-length, GC, N, and mean quality checks for locally available files.
- `results/rnaseq/phase7d_salmon_preflight.tsv`: Salmon executable, transcript FASTA, `tx2gene`, index, and paired-FASTQ readiness checks.
- `results/rnaseq/phase7d_salmon_quant_summary.tsv`: per-run Salmon execution status, command, logs, and `quant.sf` path.
- `results/rnaseq/phase7d_candidate_expression_matrix.tsv`: candidate/run-level TPM and read-count summaries only when quantification is complete.
- `reports/generated/phase7d_rnaseq_quantification_report.md`: human-readable Phase 7d status report.

If local reads, Salmon, or the Salmon index are not available, the phase writes explicit `NOT_RUN` statuses. Successful quantification can support cautious retina-specific detection summaries only. It does not support differential expression, pathway activity, whole-organism aging interpretation, or longevity-mechanism claims.

## Phase 5b Optional Repeat Annotation Manifest Notes

Phase 5b records how de novo repeat annotation should be run outside the default workflow. It does not download new public data and does not execute RepeatModeler or RepeatMasker automatically.

The planning/provenance outputs are:

- `results/repeats/phase5b_repeat_annotation_preflight.tsv`
- `results/repeats/phase5b_repeat_annotation_plan.tsv`
- `data/metadata/phase5b_repeat_annotation_provenance.tsv`
- `data/metadata/phase5b_repeat_output_inventory.tsv`

The current preflight reports the Tokyo genome FASTA at `data/raw/references/SMIC_TOKYO_GENOME_2025/GCA_056099535.1_ASM5609953v1_genomic.fna.gz` as present and readable. The external RepeatModeler/RepeatMasker route produced:

- `results/repeats/phase5b/smic_tokyo-families.fa`
- `results/repeats/phase5b/smic_tokyo.repeatmasker.out.gff`
- `results/repeats/phase5b/smic_tokyo.repeatmasker.out`

The expected `.tbl` summary was not produced, and the RepeatMasker log retains a `ProcessRepeats` signal-9 warning. These outputs are therefore used for candidate-locus artifact/context review only, not genome-wide repeat percentages. After any future external rerun, rerun the Phase 5b planning command to record checksums and file sizes, then rerun Phase 5 and Phase 5c.

Statuses used by the MVP:

- `REGISTERED`: known public resource, not downloaded.
- `LOCAL_AVAILABLE`: expected local file exists.
- `MISSING_LOCAL`: expected local file is absent.
- `NOT_ASSESSED`: resource or field has not been evaluated.
- `DEMO_ONLY`: artificial test fixture.
