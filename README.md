# AmAgo
This repository contains a code used for data analysis in a scientific paper "Prokaryotic Argonaute nuclease cooperates with co-encoded RNase to acquire guide RNAs and target invader DNA".

## What is this project for?
This is the code that was used for NGS analysis published in paper ...
1. We cleaved total RNA from *E. coli* cells with a novel RNase Agap and sequenced 5'-unphosphorylated cleavage products on Illumina platform. This directory contains scripts that process 6 NGS libraries, align the reads to the reference DNA and make a sequence Logo of Agap cleavage site.
2. We analysed small DNA and RNA associated with dAmAgo in *E. coli* cells. This directory contains scripts that process 3 NGS libraries: small RNAs, small DNAs and publicly available RNA seq data.

## 1. Agap cleavage site


Script reads_processing_and_alignment.sh starts reads quality control with FastQC, then removes adaptors with cutadapt, makes quality control once again, and then calculates the alignment.

The script has two arguments:
 
  `-p` is the number of threads to use for calculations (number of cores on the machine by default);
 
  `-d` is the path to the working directory, that should contain a single fastq.gz file and a directory that contains two FASTA files: genome.fa (the header should be >genome) and plasmid.fa (the header should be >plasmid).
 
Script cleavage_site_logo.py contains processes alignment data from Agap-treated RNA and control samples to create sequence logos for Agap cleavage sites. The script loads data, processes it, identifies cleavage sites, and generates sequence logos for visualization.

## 2. NAs associated with dAmAgo

Scripts dna_preprocessing.sh and rna_preprocessing.sh make reads quality control check with FastQC, remove adaptors with cutadapt, repet the quality control for the trimmed reads.

The script has two arguments:
 
  `-p` is the number of threads to use for calculations (number of cores on the machine by default);
 
  `-d` is the path to the working directory, that should contain a single fastq.gz file.

Scripts rna_pipeline.sh, rna_seq_pipeline.sh and dna_pipeline.sh align the reads on the plasmid and then map the analigned once to the genome. The scripts then calculate coverage of genes and of each chromosome/plasmid position.

The script has two arguments:
 
  `-p` is the number of threads to use for calculations (number of cores on the machine by default);
 
  `-d` is the path to the working directory, that should contain a single trimmed.fastq.gz file, an annotation GTF file and a directory that contains two FASTA files: genome.fa (the header should be >genome) and plasmid.fa (the header should be >plasmid). In the case of rna_seq.sh there is no need in plasmid.fa file.

Scripts logo.py, plasmid_cov.py, whole_genome_cov.py, genes_cov.py, rna_dna_cov.py lengths.py and rna_seq.py produce corresponding plots.

## Requirements
This script utilizes some commonly used programs for data analysis and NGS analysis:
- FastQC https://github.com/s-andrews/FastQC
- cutadapt https://github.com/marcelm/cutadapt
- bowtie https://github.com/BenLangmead/bowtie
- samtools https://github.com/samtools/samtools
- pandas, numpy and matplotlib libraries for Python
- pysam library for Python
- logomaker library for Python

 ## Where to find the data that was processed with this code?
 The raw sequencing reads are deposited in SRA in the BioProjects PRJNA1023518 and PRJNA1023528.
