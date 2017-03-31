# mRNA-Seq pipeline 
This pipeline performs trimming, mapping of mRNA-Seq Illumina paired-end reads on the human genome. It computes raw counts for genes and repeats and output bigwig coverage files for each sample.

# How to use it

## Install all softwares and packages needed with the [Conda package manager](https://conda.io/docs/using/envs.html)

### Create a virtual environment named "myrnaseqworkflow" from the `environment.yaml` file
conda env create --name myrnaseqworkflow --file environment.yaml

### Activate this virtual environment
source activate myrnaseqworkflow

### execute the workflow (here dry run with the -n option)
snakemake -np
 
# Main outputs:
*  FastQC quality checks (before and after trimming)
*  Mapping statistics ("logs") for each sample
*  Counts for both genes and repeats   
*  bigWig coverage file (indexed binary format of bedGraph files)

# config.yaml
This file is used to customize the analysis. 

## Directories:
*  fastqdir: this is the directory that contains the mRNA-Seq fastq files
*  workdir: this is the directory that will contain all intermediate files (can be cleaned up after analysis)
*  resultdir: this is where you want to store all desired outputs.

## Genomic references and annotations
refseqs:
  genome2bit: '../../data/02.refs/hg19.2bit'
  genomefasta: '../../data/02.refs/hg19.genome.fa'
  chromSizes: 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes'
  repeatsfasta: '../../data/02.refs/repeats.parsed.fasta'
annotations:
  gtf4genes: '../../data/02.refs/knownGenes.gtf'
  gtf4repeats: '../../data/02.refs/repeats.gtf'
  refseqFromSubread: '../../data/02.refs/hg19_RefSeq_exon.txt'

## samples
Provide a list of samples
samples:
  mDA_D0_1:
    forward: huESC_fj74_R1.fastq
    reverse: huESC_fj74_R2.fastq
  mDA_D0_2:
    forward: fj95_TGACCA_R1.fastq
    reverse: fj95_TGACCA_R2.fastq
  etc.   

## Parameters for softwares used
Trimmomatic (trimming and quality check)
STAR (mapping)
FastQC (reports of quality)

