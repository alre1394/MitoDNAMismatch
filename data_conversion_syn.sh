#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J EX_info_conversion_alre1394
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH -o "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.out"
#SBATCH -e "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.err"

module load bioinfo-tools
module load python/3.9.5

#DIRECTORIES
OUTDIR=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/annotation/extraction/syn_nonsyn
SAMPLES="$OUTDIR/samples.txt"
GENO_SYN_TSV="$OUTDIR/genotypes_for_load_syn.tsv"
LOAD_TABLE_SYN="$OUTDIR/syn_per_sample.tsv"

echo "File checks:"
wc -l "$SAMPLES"
wc -l "$GENO_SYN_TSV"
ls -lh "$GENO_SYN_TSV"
echo

echo "File checks:"
wc -l "$SAMPLES"
wc -l "$GENO_SYN_TSV"
ls -lh "$GENO_SYN_TSV"
echo

export OUTDIR

python3 << 'EOF'
import os
import sys

outdir = os.environ["OUTDIR"]
samples_file = os.path.join(outdir, "samples.txt")
syn_geno_file = os.path.join(outdir, "genotypes_for_load_syn.tsv")
syn_load_out = os.path.join(outdir, "syn_per_sample.tsv")

# Load samples
with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]

n = len(samples)
loads = [0] * n

sys.stderr.write(f"Loaded {n} samples.\n")
sys.stderr.write(f"Reading genotypes from: {syn_geno_file}\n")

line_count = 0
with open(syn_geno_file) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue

        parts = line.split("\t")
        gts = parts[2:]

        if len(gts) != n:
            sys.stderr.write(
                f"WARNING: mismatch #GT={len(gts)} vs samples={n} "
                f"at variant {line_count+1}\n")
        for i, gt in enumerate(gts):
            if i >= n:
                break
            if gt in ("./.", ".|.", ".", ""):
                continue

            alleles = gt.replace("|", "/").split("/")
            alt_count = sum(1 for a in alleles if a not in ("0", "."))
            loads[i] += alt_count

        line_count += 1
        if line_count % 10000 == 0:
            sys.stderr.write(f"Processed {line_count} variants...\n")

sys.stderr.write(f"Finished processing {line_count} variants.\n")

with open(syn_load_out, "w") as out:
    out.write("Sample\tSynonymousCounts\n")
    for s, l in zip(samples, loads):
        out.write(f"{s}\t{l}\n")

sys.stderr.write(f"Output written to: {syn_load_out}\n")
EOF
