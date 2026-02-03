#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -J IMPACT_TOTAL_EX_alre1394
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH -o "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.out"
#SBATCH -e "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.err"

module load bioinfo-tools
module load bcftools/1.20

BASE=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/annotation
AUTO_ANN="$BASE/auto_annotation.vcf.gz"
CHRX_ANN="$BASE/chrx_annotation.vcf.gz"
SAMPLES="$BASE/samples.txt"
OUTDIR="$BASE/effect_counts/"
mkdir -p "$OUTDIR"

declare -A FILTERS=(
  [HIGH]='INFO/ANN~"\|HIGH\|"'
  [MODERATE]='INFO/ANN~"\|MODERATE\|" && INFO/ANN!~"\|HIGH\|"'
  [LOW]='INFO/ANN~"\|LOW\|" && INFO/ANN!~"\|HIGH\|" && INFO/ANN!~"\|MODERATE\|"'
  [MODIFIER]='INFO/ANN~"\|MODIFIER\|" && INFO/ANN!~"\|HIGH\|" && INFO/ANN!~"\|MODERATE\|" && INFO/ANN!~"\|LOW\|"'
)

for impact in LOW MODERATE HIGH MODIFIER
do
        #AUTOSOME FILTERING
        echo "Filtering autosomes + XPAR for ${impact} impact SNPs..."
        bcftools view \
		-i "${FILTERS[$impact]}" \
		"$AUTO_ANN" \
		-Oz -o "$OUTDIR/autosomes_${impact}.vcf.gz"
	bcftools index -f "$OUTDIR/autosomes_${impact}.vcf.gz"
	
	echo "Number of autosome+XPAR ${impact} impact variants:"
        bcftools view -H "$OUTDIR/autosomes_${impact}.vcf.gz" | wc -l
        echo

        #CHRX NONPAR FILTERING
        echo "Filtering chrx NONPAR for ${impact} impact SNPs..."
        bcftools view \
		-i "${FILTERS[$impact]}" \
		"$CHRX_ANN" \
		-Oz -o "$OUTDIR/chrx_${impact}.vcf.gz"
        bcftools index -f "$OUTDIR/chrx_${impact}.vcf.gz"

        echo "Number of chrx NONPAR ${impact} impact variants:"
        bcftools view -H "$OUTDIR/chrx_${impact}.vcf.gz" | wc -l
        echo

        #MERGING impact SNPs FROM AUTO AND CHRX NONPAR
        echo "Concatenating ${impact} impact variants (autosomes + chrX NONPAR)..."
        bcftools concat -a \
            "$OUTDIR/autosomes_${impact}.vcf.gz" \
            "$OUTDIR/chrx_${impact}.vcf.gz" \
            -Oz -o "$OUTDIR/merged_${impact}.vcf.gz"
        bcftools index -f "$OUTDIR/merged_${impact}.vcf.gz"

        echo "Merged synonymous VCF: $OUTDIR/merged_${impact}.vcf.gz"
        echo "Number of merged ${impact} impact variants:"
        bcftools view -H "$OUTDIR/merged_${impact}.vcf.gz" | wc -l
        echo

        #OUTPUT TABLES FOR EVERY IMPACT SNPs (LOW, MEDIUM, HIGH, MODIFIER)
        echo "Extracting sample IDs..."
        bcftools query -l "$OUTDIR/merged_${impact}.vcf.gz" > "$SAMPLES"
        NSAMP=$(wc -l < "$SAMPLES")
        echo "Number of samples: $NSAMP"

        #CREATING GENOTYPE TABLES
        echo "Creating genotype table for ${impact} SNPs"
        bcftools query \
                -f '%CHROM\t%POS[\t%GT]\n' \
                "$OUTDIR/merged_${impact}.vcf.gz" \
                > "$OUTDIR/load_tables/genotypes_for_load_${impact}.tsv"

        echo "Number of lines (variants) in genotype table:"
        wc -l "$OUTDIR/load_tables/genotypes_for_load_${impact}.tsv"
        echo
done

echo "bcftools-only job finished at: $(date)"
