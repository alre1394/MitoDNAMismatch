#!/bin/bash -l

#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH -J NG_SNP_annotation_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH -o "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.out"
#SBATCH -e "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.err"

#IMPORTING MODULES
module load bioinfo-tools
module load snpEff/5.2
module load bcftools/1.20

#PATHS
BASE=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project
AUTO_VCF="$BASE/lukas/vcf_filtering/final_dog_mito_vcf/dog_mitonuclear.noCoyotes.vcf.gz"
CHRX_VCF="$BASE/lukas/vcf_filtering_chrx/final_dog_mito_vcf_chrx/chrx_dog_mitonuclear.noCoyotes.vcf.gz"
DIR_DB="$BASE/SNPEff/canFam4_Y.NCBIAnnotation_SNPeffV3"
CONFIG_DB="$DIR_DB/snpEff.config"
ID_DB=canFam4_y
OUTPUT="$BASE/alre/annotation"

#CHANGING DIRECTORY TO DB_DIR
cd "$DIR_DB"

#SNP ANNOTATION
#AUTOSOME GENES
snpEff -c "$CONFIG_DB" -dataDir "$DIR_DB" -v "$ID_DB" -csvStats "$OUTPUT" "$AUTO_VCF" \
	|  bcftools view -Oz --threads 8 -o "$OUTPUT/auto_annotation.vcf.gz"

bcftools index --threads 8 -f "$OUTPUT/auto_annotation.vcf.gz"

#CHRX GENES
snpEff -c "$CONFIG_DB" -dataDir "$DIR_DB" -v "$ID_DB" -csvStats "$OUTPUT" "$CHRX_VCF" \
	|  bcftools view -Oz --threads 8 -o "$OUTPUT/chrx_annotation.vcf.gz"

bcftools index --threads 8 -f "$OUTPUT/chrx_annotation.vcf.gz"

# ---- quick sanity checks ----
echo "Main VCF sample count:"
bcftools query -l "$OUTPUT/auto_annotation.vcf.gz" | wc -l

echo "chrX VCF sample count:"
bcftools query -l "$OUTPUT/chrx_annotation.vcf.gz" | wc -l

echo "Check ANN header (autosomes):"
bcftools view "$OUTPUT/auto_annotation.vcf.gz" | head -40 | grep -i ANN || echo "No ANN header found?"

echo "First annotated variant (autosomes):"
bcftools view -H "$OUTPUT/auto_annotation.vcf.gz" | head -1
