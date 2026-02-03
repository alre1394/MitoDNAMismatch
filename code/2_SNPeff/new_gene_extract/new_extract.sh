#!/bin/bash
#SBATCH -A uppmax2025-2-119
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -J EXTRACT_SNP_SAMP_alre1394
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH -o "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.out"
#SBATCH -e "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/out/%j.err"

module load bioinfo-tools
module load bcftools/1.20

BASE=/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/alre/data/annotation/new_gene_vcf
VCF_FILE="$BASE/subset.merge.vcf.gz"  # Your VCF file
SAMPLES_FILE="$BASE/samples.txt"  # File containing sample IDs

bcftools view -i 'INFO/AF<0.05 || "INFO/ANN~\"MODERATE\" || "INFO/ANN~\"HIGH\"' \
	-Oz -o filtered.gene_load.vcf.gz subset.merge.vcf.gz
(
  # Print the header
  echo -e "Sample_ID\tGene_Load"
  
  # Compute weighted gene load
  bcftools query -f '[%SAMPLE\t%GT\n]' filtered.gene_load.vcf.gz \
  | awk 'BEGIN { FS="\t" } 
    # Count hets and homs for the sample
    $2=="0/1" || $2=="1/0" || $2=="0|1" || $2=="1|0" { load[$1]+=1 }
    $2=="1/1" || $2=="1|1" { load[$1]+=2 }
    # Print final result
    END { 
      for (s in load) 
        print s "\t" load[s]
    }' > gene_load_raw.tsv
)

awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {a[$1]=0; next}
     {a[$1]=$2}
     END {
       for (s in a) 
         if (a[s] != "" && s != "") print s, a[s]
     }' \
     samples.txt gene_load_raw.tsv \
| sort > gene_load_per_sample.tsv
