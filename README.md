# Mito-nuclear mismatch
The repository stores the main scripts utilized during the project conducted in the Applied Bioinformatics course. An exploratory analysis was conducted to study the inverse correlation between dog lifespan and body weight, particularly in relation to mitochondrial dysfunction. The genetic load was quantified in mitochondrial genes located in the nuclear genome (mito-nuclear genes), based on SNP annotation and predicted impact using SnpEff. The project leveraged publicly available data from the Dog10K project [1], including 2,000 whole-genome sequenced samples.

## Overall aim:
Study the association between dog body weight and genetic load in mito-nuclear genes across different dog breeds.

## Specific aims: 
  - Filter nuclear Variant Call Format (VCF) files containing genetic variation data, in the form of single nucleotide polymorphisms (SNPs).
  - Annotate genetic mutations based on their predicted effect type and putative impact using an existing annotation tool (SnpEff).
  - Quantify genetic load per dog individual by extracting and counting "highly deleterious" SNPs (HIGH/MEDIUM impact) per dog sample.
  - Conduct an exploratory data analysis of genetic load estimations using dog body weight, as a lifespan proxy.

### Code specifics:
The code intends to perform the following tasks:
- Filter large VCF files tracking SNPs in the nuclear genome and predict their effect on gene features using SnpEff.
- Count the amount of SNPs (homozygote and heterozygotes) per dog breed in order to quantify the genetic load on a per-sample basis.
- Explore mutation data, compare different dog weight classes across weight classes, and perform a correlation analysis between dog body weight and genetic load metrics.

## Input:
- Nuclear VCF files storing SNPs per individual dog sample, relative to a dog reference genome (German Shepherd), including information such as specific genomic location (chromosome, start and end position) and allelic variants.

## Output:
- Filtered VCF files obtained with bcftools. The file excludes non-mitochondrial genes and any coyote sample, due their high variation.
- Annotated VCF file obtained with SnpEff. The tool appends an annotation column with the predicted effect type, putative impact, and allelic variant of each SNP.
- Genotype-level tsv files, containing SNP counts per dog sample.

## Code list and description:
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
    - **Data Extraction - the annotated VCF file is filtered by the respective annotation label/tag (HIGH/MEDIUM impact, synonymous or nonsynonymous effect type)**:
      - [`general_impact_extract.sh`](./code/2_SNPeff/general_impact_extract.sh)
      - [`syn_extract.sh`](./code/2_SNPeff/syn_extract.sh)
      - [`nonsyn_extract.sh`](./code/2_SNPeff/nonsyn_extract.sh)
    - **Data Conversion - creation of genotype-level tables, it counts the number of SNPs per dog sample in the filtered VCF files**:
      - [`data_conversion_impact.sh`](./code/2_SNPeff/data_conversion_impact.sh)
      - [`syn_extract.sh`](./code/2_SNPeff/syn_extract.sh)
      - [`data_conversion_nonsyn.sh`](./code/2_SNPeff/data_conversion_nonsyn.sh)
- **Data Analysis and Visualization study of genetic variation and genetic load metrics:** [`genetic_load.sh`](./code/2_SNPeff/genetic_load.sh)

## References
1.- Meadows, J.R.S., Kidd, J.M., Wang, G., et al. (2023) ‘Genome sequencing of 2000 canids by the Dog10K consortium advances the understanding of demography, genome function and architecture’, Genome Biology, 24, p. 187. https://doi.org/10.1186/s13059-023-03023-7.
