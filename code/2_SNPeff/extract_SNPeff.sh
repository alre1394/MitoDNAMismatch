#!/usr/bin/env bash

#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH -J NG_ANN_extraction_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH -o "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.out"
#SBATCH -e "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.err"

#IMPORTING MODULES
module load bioinfo-tools
module load bcftools/1.20
BASE=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/annotation
AUTO_ANN="$BASE/auto_annotation.vcf.gz"
CHRX_ANN="$BASE/chrx_annotation.vcf.gz"
OUTDIR="$BASE/extraction"

for VCF in "$AUTO_ANN" "$CHRX_ANN"; do

    PREFIX=$(basename "$VCF" .vcf.gz)
    echo "Processing: $VCF"

    # Step 1 flatten per-sample lines
    bcftools query \
        -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%INFO/ANN\n]' \
        "$VCF" \
    | awk '$6 != "0/0" && $6 != "0|0"' \
    > "$OUTDIR/${PREFIX}.all_samples_variants.tsv"

    # Step 2 extract sample + impact
    awk -F'\t' '
    {
        split($7, ann, "|");
        print $1"\t"ann[3];
    }' "$BASE/${PREFIX}.all_samples_variants.tsv" \
    > "$OUTDIR/${PREFIX}.sample_impact.tsv"

    # Step 3 count impacts per sample
    awk '
    {
        s=$1; imp=$2;
        if (imp=="HIGH") high[s]++;
        else if (imp=="MODERATE") moderate[s]++;
        else if (imp=="LOW") low[s]++;
    }
    END {
        print "Sample\tHIGH\tMODERATE\tLOW";
        for (s in high)
            printf "%s\t%d\t%d\t%d\n", s, high[s], moderate[s], low[s];
    }' \
    "$OUTDIR/${PREFIX}.sample_impact.tsv" \
    > "$OUTDIR/${PREFIX}.impact_summary.tsv"

done
