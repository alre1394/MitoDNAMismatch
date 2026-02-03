#!/usr/bin/env bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J TOTAL_SAMP_COUNT_alre1394
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH -o "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.out"
#SBATCH -e "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.err"

module load bioinfo-tools
module load bcftools/1.20

BASE="/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/annotation"
VCF="$BASE/merged_mt.vcf.gz"
SAMP="$BASE/samples.txt"
OUT="$BASE/count_samp.tsv"

#File header
echo -e "sample_id\tsnp_count" > "$OUT"

# Count SNPs per sample (non-reference genotypes only) considering all genotypes per samples
bcftools view -v snps "$VCF" \
| bcftools query -f '[%SAMPLE\t%GT\n]' \
| awk '
$2!="0/0" && $2!="0|0" && $2!="./." && $2!=".|." {
    count[$1]++
}
END {
    for (s in count)
        print s, count[s]
}' OFS="\t" \
| sort >> "$OUT"
