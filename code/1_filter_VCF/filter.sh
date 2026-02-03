#!/bin/bash -l

#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J NG_VCF_filter_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH -o "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.out"
#SBATCH -e "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.err"

#IMPORTING MODULES
module load bioinfo-tools
module load bcftools/1.9

#INPUT DATA
BASE=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project
VCF="$BASE/data/ng_VCFs"
AUTO_VCF="$VCF/AutoAndXPAR.SNPs.vcf.vqsr99.gz"
CHRX_VCF="$VCF/chrX.NONPAR.SNPs.vcf.vqsr99.gz"
BED_NORM="$BASE/alre/filter_ng/filtering_vcf/mito_nuclear_genes_chr_tab.sorted.bed"
COYOTE_IDS="$BASE/alre/filter_ng/filtering_vcf/coyote.txt"

#BCFtools FILTERING AUTOSOMES
bcftools view \
	--regions-overlap pos \
	-S ^"$COYOTE_IDS" \
	-R "$BED_NORM" -Oz "$AUTO_VCF" -o auto_ng_filter.vcf.gz
bcftools index auto_ng_filter.vcf.gz

#BCFtools FILTERING CHRX
bcftools view -S ^"$COYOTE_IDS" -R "$BED_NORM" -Oz "$CHRX_VCF" -o chrx_ng_filter.vcf.gz
bcftools index chrx_ng_filter.vcf.gz
