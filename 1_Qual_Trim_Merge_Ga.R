# Combinatorial library quality filtering, trimming, and read merging pipeline

# Author: Steven Chen
# Email: stevenkuanyu.chen@gmail.com
# Version History
#	2023-07-12: v3

# =============================================================================
# Preamble
# =============================================================================

# Summary:
# Cutadapt is used for trimming. Sequences used to inform trimming (primer inputs)
# were rationalized a priori by looking at quality control checks of the reads in FASTQC.
# Parts of the DADA2 pipeline (https://benjjneb.github.io/dada2/tutorial.html) is used
# for light trimming and quality filtering. PANDAseq will then merge the read pairs 
# (using either its default simple_bayesian or UPARSE/USEARCH algorithms).

# Programs and packages summary:
# R v4.2.2 (https://www.r-project.org)
# DADA2 v1.26.0 (https://www.nature.com/articles/nmeth.3869)
# Cutadapt v4.2 (https://journal.embnet.org/index.php/embnetjournal/article/view/200)
# Pandaseq v2.11 (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-31)

# To run this script, make sure to:
# 1. Install Anaconda (https://www.anaconda.com) or Miniconda (https://docs.conda.io/en/latest/miniconda.html).
# 2. Create a conda environment (https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
# and install dada2 package by executing `conda install -c bioconda bioconductor-dada2`in the linux terminal.
# 3. In the same conda environment, install cutadapt and pandaseq and make sure that the script below is modified to
# point to their installed locations.
# 4. Create a directory containing this R script and the relevant FASTQ files pertaining to a single variant library
# (Ga OR H3; A separate directory, R script, and set of demultiplexed FASTQ files should be used for the other library)
# 5. Set the appropriate work directories in the script below
# 6. Set appropriate primer sequences for cutadapt (variant library-dependent). For Ga use "CTATGAAGTTCGTATTGAGT" for
# FWD and "TAAGAAATTCGCGCGGCCGCTTATTTA" for REV. For H3 use "GGATGCAATTTGGAGGGCTT" for FWD and 
# "AGTTGCTCATGGGCTTACACACCACCA" for REV.
# 7. To execute this script, navigate to the directory containing this script in the terminal and
# execute `R --no-save < [NAME OF .R script] > [NAME OF LOG FILE].log`

# This script might take a couple hours to complete depending on the computer.


# =============================================================================
# Quality filtering and trimming
# =============================================================================

# Set work directory
setwd("/media/scratch/stevenscratch/Comb_lib_seq_data/2023-07-12_Ga+H3_trimQC_UPARSE/Ga")

# print date and time of run for .log file
Sys.time()

# Activate libraries and packages
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(stringr)
packageVersion("stringr")

# Set path
path <- "/media/scratch/stevenscratch/Comb_lib_seq_data/2023-07-12_Ga+H3_trimQC_UPARSE/Ga"

list.files(path) # to check that the packages have loaded successfully

# Perform some string manipulation to get matched list of the forward and reverse fastq files
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))

# Identify primers for cutadapt trimming (can also be used to remove low-quality base positions)
FWD <- "CTATGAAGTTCGTATTGAGT"  ## Gaf... CHANGE ME FOR DIFFERENT PRIMERS
REV <- "TAAGAAATTCGCGCGGCCGCTTATTTA"  ## Gar... CHANGE ME

# Check primer orientations
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Filter for ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # make sure all files are .fastq.gz or will throw error

# As a sanity check, we will count the number of times the primers appear in the forward and reverse reads:
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# =============================================================================
# Cutadapt to removing primers and flanking low-quality bases
# =============================================================================

# Initiating Cutadapt
cutadapt <- "/home/steven/anaconda3/envs/amplicon/bin/cutadapt" # CHANGE ME to the cutadapt location on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

# Executing Cutadapt: create output filenames for Cutadapt-ed files, define the paramaters, and execute cutadapt to remove the primers
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC, "-m", 1) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC, "-m", 1) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# As a sanity check, we will count the presence of primers in the first cutadapt-ed sample:
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Read in the names of the cutadapt-ed FASTQ files and applying some string manipulation to get the matched lists of forward and reverse fastq files
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_L")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Filter and trim

# Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# Setting parameters and conducting filtering and trimming
# parameter for Ga
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(1, 1),
    truncQ = 11, minLen = 64, maxLen = 64, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)

# =============================================================================
# Pandaseq merging
# =============================================================================

pandaseq <- "/home/steven/anaconda3/envs/amplicon/bin/pandaseq" # CHANGE ME to the pandaseq location on your machine
system2(pandaseq, args = "-h") # Run shell commands from R

cut.filt.path <- "/media/scratch/stevenscratch/Comb_lib_seq_data/2023-07-12_Ga+H3_trimQC_UPARSE/Ga/cutadapt/filtered" # CHANGE ME to the path where you want to create folder for merged reads

# Read in the names of the Cutadapt-ed + filtered files
cut.filtFs <- sort(list.files(cut.filt.path, pattern = "_R2_001.fastq.gz", full.names = TRUE))
cut.filtRs <- sort(list.files(cut.filt.path, pattern = "_R1_001.fastq.gz", full.names = TRUE))

# Execute Pandaseq
path.merge <- file.path(cut.filt.path, "pandamerge")
if(!dir.exists(path.merge)) (dir.create(path.merge))
filtFs.merge <- file.path(path.merge, basename(cut.filtFs))
# filtRs.merge <- file.path(path.merge, basename(cut.filtRs))

for(i in seq_along(cut.filtFs)) {
    system2(pandaseq, args = c("-A", "uparse", "-F", "-f", cut.filtFs[i],
                               "-r", cut.filtRs[i], "-d",
                               "rbfkms", "-o", 63, "-O", 65,
                               "-W", filtFs.merge[i]))
}

# Pandaseq outputs .bz2 files but the basename used to generate the
# merged files follows the .fastq.gz so this is to rename it to .bz2
pandamerge.path <- "/media/scratch/stevenscratch/Comb_lib_seq_data/2023-07-12_Ga+H3_trimQC_UPARSE/Ga/cutadapt/filtered/pandamerge" # CHANGE ME to directory containing merged files

setwd(pandamerge.path) # set work directory, required to rename files

# Execute file renaming
file.rename((list.files(pandamerge.path, pattern = ".gz")), 
            (str_replace(list.files(pandamerge.path, pattern = ".gz"), pattern = ".gz", ".bz2")))
mergedF.bz2 <- list.files(pandamerge.path, pattern = ".bz2")
# Unzip the bzip2 files
bzip2 <- "/home/steven/anaconda3/envs/amplicon/bin/bzip2" # Path to bzip2 files
system2(bzip2, "--version")
# Execute bzip2
system2(bzip2, args = c("-d", "*.fastq.bz2"))

# =============================================================================
# Additional useful commands
# =============================================================================

# The merged files can then be concatenated into a single file for processing
# by using the `cat` command (available in Linux/Unix). 
# Example: `cat *.fastq > Ga_merged.fastq`

# Convert from concatenated merged.fastq file to .fasta file for analysis.
# Example: `sed -n '1~4s/^@/>/p;2~4p' INFILE.fastq > OUTFILE.fasta`

# Quantify the number of unique sequences using Mothur.
# Example: mothur > `unique.seqs(fasta = INFILE.fasta)`
