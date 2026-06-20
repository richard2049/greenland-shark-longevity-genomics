# Candidate Domain Validation

This phase checks representative candidate proteins against Pfam-A with HMMER. It is a protein-model validation step, not a biological conclusion.

On the current native Windows environment, the workflow uses `pyhmmer 0.12.1` because command-line `hmmscan` is not available through the current `win-64` conda channels. This is still a HMMER/Pfam profile-HMM scan. Linux, WSL, or container users can switch `domain_validation.backend` to `hmmscan` if a command-line HMMER installation is available.

## Inputs

- Candidate FASTA: `data/interim/candidate_validation/representative_candidate_proteins.faa`
- Domain targets: `results/validation/domain_check_targets.tsv`
- Pfam database: `data/raw/references/PFAM/Pfam-A.hmm`

Pfam-A can be obtained from the Pfam current-release FTP area. The current local installation is:

- Release: `38.1`
- Retrieval date: `2026-05-29`
- URL: `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`
- `Pfam-A.hmm.gz` SHA-256: `D3D30C8E6801BFEDECF783408ECC98916F8F1DDA8974C6E51036FCBDD765F591`
- `Pfam-A.hmm` SHA-256: `A78A42D6FAF265B6BFCA59E8F062D06FAE6083CE2C6E335D7B381F20B82B7903`

The expected compressed file is commonly distributed as:

```text
https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
```

If using command-line `hmmscan`, run:

```bash
hmmpress data/raw/references/PFAM/Pfam-A.hmm
```

The current `pyhmmer` backend can scan the plain HMM file directly and does not require `hmmpress` indexes.

## Workflow

Run:

```bash
snakemake --snakefile workflow/Snakefile --cores 4 --config workflow_mode=domain_validation
```

This mode first rebuilds the real orthology-postprocess targets. With the current `pyhmmer` backend it then performs a hmmscan-compatible profile-HMM scan and writes HMMER-style domain-table output. The equivalent command-line route is:

```bash
hmmscan --cpu 4 --cut_ga --domtblout results/domains/representative_candidate_proteins.pfam.domtblout --noali data/raw/references/PFAM/Pfam-A.hmm data/interim/candidate_validation/representative_candidate_proteins.faa
```

The configured command is stored in `config/config.yaml`.

## Outputs

- `results/validation/hmmer_pfam_preflight.tsv`
- `results/domains/representative_candidate_proteins.pfam.domtblout`
- `results/domains/representative_candidate_proteins.pfam.txt`
- `results/validation/candidate_domain_integrity.tsv`

## Classification

The parser classifies each representative candidate protein as:

- `DOMAIN_SUPPORTED`: at least one Pfam hit passes the configured independent e-value and full-domain HMM coverage thresholds.
- `PARTIAL_DOMAIN`: a significant Pfam hit is present, but full-domain coverage is incomplete.
- `NO_EXPECTED_DOMAIN_DETECTED`: no qualifying Pfam hit is detected under the configured thresholds.
- `NOT_ASSESSED`: the representative protein or HMMER/Pfam run is unavailable.

These are domain-evidence labels only. They do not prove functional activity, gene duplication validity, inactivation, absence, or longevity relevance.

## InterProScan

InterProScan should be treated as a later confirmation and annotation expansion layer. It is broader than HMMER/Pfam because it integrates multiple InterPro member databases, but that breadth is not required for the first candidate-domain integrity check.

## Downstream Annotation Rescue

The next workflow mode is `annotation_rescue`, which targets unresolved high-priority genes left as annotation uncertainty after exact-symbol orthology mapping. It uses reciprocal protein similarity to prioritize focal protein models and writes the Phase 4 rescue tables under `results/rescue/`. These outputs remain protein-level rescue candidates and require genome-level alignment before annotation, duplication, or loss language is justified.

The follow-up mode is `genome_validation`. It prepares miniprot queries and runs splice-aware protein-to-genome alignment only when the target genome FASTA and miniprot executable are available. On the current laptop this runs through `tools/miniprot-docker.cmd` in `candidate_scaffolds` mode. The output remains candidate genome-alignment evidence, not a final biological claim.
