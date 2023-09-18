# Comb_lib_NGS_analysis
[![DOI](https://zenodo.org/badge/661528862.svg)](https://zenodo.org/badge/latestdoi/661528862)

## Part 1: Pipeline for quality filtering, trimming, and merging
### Summary:
Cutadapt is used for trimming. Sequences used to inform trimming (primer inputs)
were rationalized a priori by looking at quality control checks of the reads in FASTQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
Parts of the DADA2 pipeline (https://benjjneb.github.io/dada2/tutorial.html) is used
for light trimming and quality filtering. PANDAseq will then merge the read pairs 
(using either its default simple_bayesian or UPARSE/USEARCH algorithms).
### Required programs and packages summary:
- R v4.2.2 (https://www.r-project.org)
- DADA2 v1.26.0 (https://www.nature.com/articles/nmeth.3869)
- Cutadapt v4.2 (https://journal.embnet.org/index.php/embnetjournal/article/view/200)
- Pandaseq v2.11 (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-31)
### Instructions for running this script:
1. Install Anaconda (https://www.anaconda.com) or Miniconda (https://docs.conda.io/en/latest/miniconda.html).
2. Create a conda environment (https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) and install dada2 package by executing `conda install -c bioconda bioconductor-dada2`in the linux terminal.
3. Activate the conda environment, install Cutadapt and Pandaseq and make sure that the script below is modified to point to their installed locations.
4. Create a directory containing this R script and the relevant FASTQ files pertaining to a single variant library (Ga OR H3; A separate directory, R script, and set of demultiplexed FASTQ files should be used for the other library).
5. Set the appropriate work directories in the script below.
6. Set appropriate primer sequences for Cutadapt (variant library-dependent). For Ga use "CTATGAAGTTCGTATTGAGT" for FWD and "TAAGAAATTCGCGCGGCCGCTTATTTA" for REV. For H3 use "GGATGCAATTTGGAGGGCTT" for FWD and "AGTTGCTCATGGGCTTACACACCACCA" for REV.
7. To execute this script, navigate to the directory containing this script in the terminal and execute `R --no-save < [NAME OF .R script] > [NAME OF LOG FILE].log`.
8. Once the run is complete, navigate into the output directory (cutadapt > filtered > pandamerge) to find the output files. These files can then be concatenated into a single file and converted into fasta format for further downstream processing (see bottom of script for commands).

Note: This script might take a couple hours to complete depending on the computer.

## Part 2: Data conversion for further processing
### Summary:
FASTA files corresponding to each variant library (Ga and H3) are converted
into csv files to enable easy separation of bases at non-conserved (mutated)
and conserved base positions. After separation, the mutated and non-mutated
bases are written into parquet file format (uses RAPIDS and/or DASK on NVIDIA
GPUs) to enable much faster read and write speeds compared to csv format.
The mutated bases are then concatenated into a single file.
### Required programs and packages summary:
- Python3.9 (https://www.python.org/downloads/)
- Rapids22.12.00 (https://docs.rapids.ai/install#conda)
- Dask-cudf22.12.01

## Part 3: Coverage analysis
### Summary:
Here, we determine the number of unique dna and protein sequences in each
library given a range of abundance filters.
### Required programs and packages summary:
Same as Part 2.


## Part 4: Variant abundance analysis
### Summary:
This script takes pre-processed sequence files (parquet files) containing
regions designated for mutagenesis (mutated) and surrounding bases
(non-mutated). To determine the proportion of false positives DNA variant
calls given all reads, the mutated regions (12 bases) were compared to 12
randomly selected bases from the non-mutated regions for each library. The
resulting distribution of variant abundances for both the mutated and
non-mutated regions can then be compared to set help abundance thresholds. 
Abundance thresholds are arbitrary but computationally necessary. We 
justify our selection here by looking for a natural break in the 
distributions. Proportion of false positive variant calls given all reads 
can then be determined for user-selected abundance thresholds.
### Required programs and packages summary:
Same as Part 2.

### Relevant References
- Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581–583. https://doi.org/10.1038/nmeth.3869 (DADA2)
- Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.Journal, 17(1), 10. https://doi.org/10.14806/ej.17.1.200 (Cutadapt)
- Masella, A. P., Bartram, A. K., Truszkowski, J. M., Brown, D. G., & Neufeld, J. D. (2012). PANDAseq: Paired-end assembler for illumina sequences. BMC Bioinformatics, 13(1). https://doi.org/10.1186/1471-2105-13-31 (PANDAseq)
- Edgar, R. C. (2013). UPARSE: highly accurate OTU sequences from microbial amplicon reads. Nature Methods 2013 10:10, 10(10), 996–998. https://doi.org/10.1038/nmeth.2604 (UPARSE/USEARCH algorithm implemented in PANDAseq)
- Rosen, M. J., Callahan, B. J., Fisher, D. S., & Holmes, S. P. (2012). Denoising PCR-amplified metagenome data. http://www.biomedcentral.com/1471-2105/13/283 (Abundance thresholding idea)
