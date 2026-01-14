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
AUTOSOME_NONSYN="$OUTDIR/autosomes_nonsyn.vcf.gz"
CHRX_NONSYN="$OUTDIR/chrx_nonsyn.vcf.gz"
MERGED_NONSYN_OUT="$OUTDIR/nonsyn_merged.vcf.gz"
NONSYN_GT_TABLE="$OUTDIR/genotypes_for_load_nonsyn.tsv"

#AUTOSOME FILTERING
echo "Filtering autosomes + XPAR for NONSYNONYMOUS SNPs..."
bcftools view \
    -i 'INFO/ANN ~ "missense_variant" || INFO/ANN ~ "stop_gained" || INFO/ANN ~ "frameshift_variant" || INFO/ANN ~ "stop_lost" || INFO/ANN ~ "stop_retained_variant"' \
    "$AUTO_ANN" \
    -Oz -o "$AUTOSOME_NONSYN"

bcftools index -f "$AUTOSOME_NONSYN"
echo "Number of autosome+XPAR nonsynonymous variants:"
bcftools view -H "$AUTOSOME_NONSYN" | wc -l
echo

#CHRX FILTERING
echo "Filtering chrX NONPAR for NONSYNONYMOUS SNPs..."
bcftools view \
    -i 'INFO/ANN ~ "missense_variant" || INFO/ANN ~ "stop_gained" || INFO/ANN ~ "frameshift_variant" || INFO/ANN ~ "stop_lost" || INFO/ANN ~ "stop_retained_variant"' \
    "$CHRX_ANN" \
    -Oz -o "$CHRX_NONSYN"

bcftools index -f "$CHRX_NONSYN"
echo "Number of chrX NONPAR nonsynonymous variants:"
bcftools view -H "$CHRX_NONSYN" | wc -l
echo

#MERGING NONSYNONYMOUS VARIANTS FROM AUTO AND CHRX NONPAR
echo "Concatenating nonsynonymous variants (autosomes + chrX NONPAR)..."
bcftools concat -a \
    "$AUTOSOME_NONSYN" \
    "$CHRX_NONSYN" \
    -Oz -o "$MERGED_NONSYN_OUT"

bcftools index -f "$MERGED_NONSYN_OUT"
echo "Merged nonsynonymous VCF: $MERGED_NONSYN_OUT"
echo "Number of merged nonsynonymous variants:"
bcftools view -H "$MERGED_NONSYN_OUT" | wc -l
echo

#OUTPUT TABLES FOR NONSYNONYMOUS VARIANTS
echo "Extracting sample IDs..."
bcftools query -l "$MERGED_NONSYN_OUT" > "$SAMPLES"
NSAMP=$(wc -l < "$SAMPLES")
echo "Number of samples: $NSAMP"

echo "Creating genotype table (this may take a bit)..."
bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$MERGED_NONSYN_OUT" > "$NONSYN_GT_TABLE"
echo "Genotype table written to: $NONSYN_GT_TABLE"
echo "Number of lines (variants) in genotype table:"
wc -l "$NONSYN_GT_TABLE"
echo "bcftools-only job finished at: $(date)"
