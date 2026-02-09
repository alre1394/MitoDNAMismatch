# Mito-nuclear mismatch
Store scripts utilized in the mito-nuclear mismatch project. An exploratory analysis was conducted to find evidence between the genetic load and life expectancy across multiple dog breeds.
Effectat the nuclear genome was quantified bases on predicted annotations

## Overall aim:
Study the association between dog body weight and genetic load in mitochondrial-related nuclear-encoded (mito-nuclear) genes.

## Specific aims: 
  - Process Variant Call Format (VCF) files containing mutation data
  - Annotate genetic mutations based on their predicted effect type and putative impact
  - Quantify genetic load per dog individual
  - Exploratory data analysis of mutation data
  - 
### Code specifics:
The code intends to perform the following tasks:
- Filter large VCF files tracking SNPs in the nuclear genome and predict their effect on gene features using SnpEff.
- Explore mutation data, compare different dog weight classesacross weight classes, and perform a correlation analysis between dog body weight and genetic load metrics.
- 
### Input:
- Nuclear VCF files storing SNPs per individual dog sample relative to a dog reference genome German sheepherd, including information such as specific genomic location (chromosome, start and end position) and allelic variants per sample.

## Output:
- Filtered VCF files obtained with bcftools. The file excludes non-mitochondrial genes.
- Annotated VCF file obtained with SnpEff. The tool appends an annotation column (8th) with the predicted effect type, putative impact, allelic variant per each SNP.

### Code list and description:
All code used for the study and belonging to each analysis is listed:

### Core Genetic Load Pipeline
- **Build target region + subset VCF:** [`vcf_filtering_translation.txt`](./code/2_SNPeff/vcf_filtering_translation.txt)  
- **Subsetting genome-wide VCF to mito-nuclear gene regions (excluding coyotes):** [`bcftools.sh`](./code/2_SNPeff/bcftools.sh)
- **Subsetting ChrX (non-PAR) VCF to mito-nuclear genes (excluding coyotes):** [`bcftools_chrx.sh`](./code/2_SNPeff/bcftools_chrx.sh)
- **Annotating mito-nuclear subset VCFs with snpEff (autosomes + XPAR and ChrX NONPAR):** [`annotate_mito_vcfs.sh`](./code/2_SNPeff/annotate_mito_vcfs.sh)
- **Extracting damaging mito-nuclear variants and generating genotype tables:** [`genetic_load.sh`](./code/2_SNPeff/genetic_load.sh) 
- **Computing per-sample genetic load from damaging variant genotypes:** [`genetic_load_python.sh`](./code/2_SNPeff/genetic_load_python.sh)

### Core Pathway-level Analysis and Data Examination:

- **Annotation of VCF file for Nuclear Autosome + Chrx-NonPar SNPs:** [`annotation_SNPeff.sh`](./code/2_SNPeff/annotation_SNPeff.sh)  
- **Pathway-level genetic load analyses:**
  - **Filtering annotated VCF files, data extraction and conversion - General Impact, Synonymous and Nonsynonymous Variant Analysis:**
    - **Data Extraction**:
      - [`general_impact_extract.sh`](./code/2_SNPeff/general_impact_extract.sh)
      - [`syn_extract.sh`](./code/2_SNPeff/syn_extract.sh)
      - [`nonsyn_extract.sh`](./code/2_SNPeff/nonsyn_extract.sh)
    - **Data Conversion**:
      - [`data_conversion_impact.sh`](./code/2_SNPeff/data_conversion_impact.sh)
      - [`syn_extract.sh`](./code/2_SNPeff/syn_extract.sh)
      - [`data_conversion_nonsyn.sh`](./code/2_SNPeff/data_conversion_nonsyn.sh)
- **Data Analysis and Visualization study of genetic variation and genetic load metrics:** [`genetic_load.sh`](./code/2_SNPeff/genetic_load.sh)

