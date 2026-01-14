#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -J ANN_info_ex_alre1394
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH -o "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.out"
#SBATCH -e "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.err"

module load bioinfo-tools
module load bcftools/1.22

BASE=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/annotation
AUTO_ANN="$BASE/auto_annotation.vcf.gz"
CHRX_ANN="$BASE/chrx_annotation.vcf.gz"
SAMPLES="$BASE/samples.txt"
OUTDIR="$BASE/extraction/syn_nonsyn"

#OUTPUT PATHS
AUTOSOME_SYN="$OUTDIR/autosomes_syn.vcf.gz"
CHRX_SYN="$OUTDIR/chrx_syn.vcf.gz"
MERGED_SYN_OUT="$OUTDIR/syn_merged.vcf.gz"
SYN_GT_TABLE="$OUTDIR/genotypes_for_load_syn.tsv"

#AUTOSOME FILTERING
echo "Filtering autosomes + XPAR for SYNONYMOUS SNPs..."
bcftools view \
    -i 'INFO/ANN ~ "synonymous_variant"' \
    "$AUTO_ANN" \
    -Oz -o "$AUTOSOME_SYN"

bcftools index -f "$AUTOSOME_SYN"
echo "Number of autosome+XPAR synonymous variants:"
bcftools view -H "$AUTOSOME_SYN" | wc -l
echo

#CHRX FILTERING
echo "Filtering chrX NONPAR for SYNONYMOUS SNPs..."
bcftools view \
    -i 'INFO/ANN ~ "synonymous_variant"' \
    "$CHRX_ANN" \
    -Oz -o "$CHRX_SYN"

bcftools index -f "$CHRX_SYN"
echo "Number of chrX NONPAR synonymous variants:"
bcftools view -H "$CHRX_SYN" | wc -l
echo

#MERGING SYNONYMOUS VARIANTS FROM AUTO AND CHRX NONPAR
echo "Concatenating synonymous variants (autosomes + chrX NONPAR)..."
bcftools concat -a \
    "$AUTOSOME_SYN" \
    "$CHRX_SYN" \
    -Oz -o "$MERGED_SYN_OUT"

bcftools index -f "$MERGED_SYN_OUT"
echo "Merged synonymous VCF: $MERGED_SYN_OUT"
echo "Number of merged synonymous variants:"
bcftools view -H "$MERGED_SYN_OUT" | wc -l
echo

#OUTPUT TABLES FOR SYNONYMOUS AND NONSYNONYMOUS VARIANTS
echo "Extracting sample IDs..."
bcftools query -l "$MERGED_SYN_OUT" > "$SAMPLES"
NSAMP=$(wc -l < "$SAMPLES")
echo "Number of samples: $NSAMP"

echo "Creating genotype table (this may take a bit)..."
bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$MERGED_SYN_OUT" > "$SYN_GT_TABLE"
echo "Genotype table written to: $SYN_GT_TABLE"
echo "Number of lines (variants) in genotype table:"
wc -l "$SYN_GT_TABLE"
echo "bcftools-only job finished at: $(date)"
