#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J NG_VCF_comparison_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH -o "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.out"
#SBATCH -e "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.err"

# VCF file paths
vcf_file1="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/filter_ng/final_vcf/chrx_ng_filter.vcf.gz"
vcf_file2="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/lukas/vcf_filtering_chrx/final_dog_mito_vcf_chrx/chrx_dog_mitonuclear.noCoyotes.vcf.gz"

# Extract Chromosome IDs and SNPs from both files (ignoring header lines starting with '#')
# Count number of SNPs for each file and list unique Chromosome IDs

# File 1 (compressed)
snps_file1=$(zcat "$vcf_file1" | grep -v '^#' | wc -l)
chromosomes_file1=$(zcat "$vcf_file1" | grep -v '^#' | awk '{print $1}' | sort | uniq)

# File 2 (compressed)
snps_file2=$(zcat "$vcf_file2" | grep -v '^#' | wc -l)
chromosomes_file2=$(zcat "$vcf_file2" | grep -v '^#' | awk '{print $1}' | sort | uniq)

# Compare number of SNPs
echo "Number of SNPs in $vcf_file1: $snps_file1"
echo "Number of SNPs in $vcf_file2: $snps_file2"

if [ "$snps_file1" -eq "$snps_file2" ]; then
    echo "Both files have the same number of SNPs."
else
    echo "The files have a different number of SNPs."
fi

# Compare Chromosome IDs
echo -e "\nChromosome IDs in $vcf_file1:"
echo "$chromosomes_file1"
echo -e "\nChromosome IDs in $vcf_file2:"
echo "$chromosomes_file2"

if diff <(echo "$chromosomes_file1") <(echo "$chromosomes_file2") > /dev/null; then
    echo "Both files have the same Chromosome IDs."
else
    echo "The Chromosome IDs differ between the two files."
fi
