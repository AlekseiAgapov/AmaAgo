# AmaAgo
This repository contains a code used for data analysis in a scientific paper ...

## What is this project for?
This is the code that was used for NGS analysis published in paper ...
We cleaved total RNA from *E. coli* cells and sequenced 5'-unphosphorylated cleavage products on Illumina platform.
This repository contains scripts that process 6 NGS libraries, align the reads to the reference DNA and make a sequence Logo of AgAP cleavage site.

## How to use it?


Script reads_processing_and_alignment.sh starts reads quality control with FastQC, then removes adaptors with cutadapt, makes quality control once again, and then calculates the alignment.

The script has two arguments:
 
  `-p` is the number of threads to use for calculations (number of cores on the machine by default);
 
  `-d` is the path to the working directory, that should include trimmed.fastq.gz file (result of reads_preprocessing.sh script) and a directory that contains three fasta files: genome.fa (the header should be >genome), plasmid.fa (the header should be >plasmid), and genome.gff3 (note that the name of the chromosome should be "genome").
 
Script cleavage_site_logo.py contains processes alignment data from AgAP-treated RNA and control samples to create sequence logos for AgAP cleavage sites. The script loads data, processes it, identifies cleavage sites, and generates sequence logos for visualization.


## Requirements
This script utilizes some commonly used programs for data analysis and NGS analysis:
- FastQC https://github.com/s-andrews/FastQC
- cutadapt https://github.com/marcelm/cutadapt
- bowtie https://github.com/BenLangmead/bowtie
- samtools https://github.com/samtools/samtools
- pandas, numpy and matplotlib libraries for Python
- logomaker library for Python

 ## Where to find the data that was processed with this code?
 The raw sequencing reads are deposited in SRA in the BioProject XXX.
