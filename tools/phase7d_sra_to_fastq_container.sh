#!/bin/sh
set -eu

run="$1"
threads="${2:-8}"
raw=/work/data/raw/rnaseq/SMIC_RETINA_PRJNA1246101_2026
workdir="$raw/_sra_work"

mkdir -p "$raw" "$workdir/prefetch" "$workdir/tmp/$run"

if [ -s "$raw/${run}_1.fastq.gz" ] && [ -s "$raw/${run}_2.fastq.gz" ]; then
  if gzip -t "$raw/${run}_1.fastq.gz" && gzip -t "$raw/${run}_2.fastq.gz"; then
    echo "[$(date)] $run final FASTQ files already exist and pass gzip integrity; skipping."
    exit 0
  fi
fi

echo "[$(date)] Starting prefetch for $run"
prefetch --max-size 100G -O "$workdir/prefetch" "$run"
sra="$workdir/prefetch/$run/$run.sra"
if [ ! -s "$sra" ]; then
  echo "Missing expected SRA file: $sra" >&2
  exit 2
fi

echo "[$(date)] Starting fasterq-dump for $run"
rm -f "$raw/${run}_1.fastq" "$raw/${run}_2.fastq" "$raw/${run}_1.fastq.gz" "$raw/${run}_2.fastq.gz"
rm -rf "$workdir/tmp/$run"
mkdir -p "$workdir/tmp/$run"
fasterq-dump --split-files --threads "$threads" --temp "$workdir/tmp/$run" --outdir "$raw" "$sra"
test -s "$raw/${run}_1.fastq"
test -s "$raw/${run}_2.fastq"

echo "[$(date)] Compressing FASTQ files for $run"
gzip -f "$raw/${run}_1.fastq"
gzip -f "$raw/${run}_2.fastq"
gzip -t "$raw/${run}_1.fastq.gz"
gzip -t "$raw/${run}_2.fastq.gz"
ls -lh "$raw/${run}_1.fastq.gz" "$raw/${run}_2.fastq.gz"
echo "[$(date)] Completed $run"
