Scripts for projects at Frank Jacobs group (University of Amsterdam, Primate Genome Evolution and Human Brain Development) 
You can clone the whole directory using: `git clone https://github.com/mgalland/neurofjacobs.git` 

## Data origin
### Annotations
The repeats.gtf, repeats.fasta and knownGenes.gtf were downloaded using the UCSC table browser (2016-09-12). Human genome __hg19/GRCh37__ was used. RepeatMasker was used to get the repeats annotation. 
### Genome sequence
The hg19.genome.fa sequence was obtained by converting the `hg19.2bit` to FASTA format using the `twoBitToFa` program (available online). 
### Sequencing data
Files for samples (Day 0) were downloaded from hgwdev.sdsc.edu (on September 9, 2016) 
*  huESC_fj74_L4_1.fastq.gz
*  huESC_fj74_L4_2.fastq.gz
*  huESC_fj74_L5_1.fastq.gz
*  huESC_fj74_L5_2.fastq.gz
*  fj95_TGACCA_L006_R1_001.fastq.gz to fj95_TGACCA_L006_R1_012.fastq.gz
*  fj95_TGACCA_L006_R2_001.fastq.gz to fj95_TGACCA_L006_R2_012.fastq.gz
*  fj96_ACAGTG_L006_R1_001.fastq.gz to fj96_ACAGTG_L006_R1_012.fastq.gz
*  fj96_ACAGTG_L006_R1_001.fastq.gz to fj96_ACAGTG_L006_R1_012.fastq.gz

The rest of the sequences were available directly from Nina Haring. 
The sample to fastq files correspondence can be found in the 'data/01.seqruns/20160912_dopamine_spheres/samples.txt' file

## Scr (source)
In this folder you will various softwares that were used for the analysis (e.g. twoBitToFa to convert hg19.2bit to fasta format)
*  To use a software, first download it or use the version in this project repository `scr/`
*  Then type in `nano ~/.bash_profile` and write `export PATH=$PATH:[/yourpathto/scr/]` to edit your $PATH variable so that it contains the path to the software. Save the `bash_profile` file (Ctrl+X) and restart a bash session and type the name of the program e.g. `twoBitToFa` to see if you can access it.

### subread
This package contains the `featureCounts` program used to compute read counts for selected level (counts at gene level). In the same folder, one can also find reference annotation (including the hg19 genome sequence) reformated specifically for featureCounts. 
>From the _subread_ documentation:
>In-built gene annotations for genomes hg19, mm10 and mm9 are included in both Bioconductor Rsubread package and SourceForge Subread package. These annotations were downloaded from NCBI RefSeq database and then adapted by merging overlapping exons from the same
>gene to form a set of disjoint exons for each gene. Genes with the same Entrez gene identiers were also merged into one gene.
>Each row in the annotation represents an exon of a gene. There are ve columns in the annotation data including Entrez gene identier (GeneID), chromosomal name (Chr),
chromosomal start position(Start), chromosomal end position (End) and strand (Strand). 

### bedGraphToBigWig
From the UCSC softwares. Used to convert bedGraph files into bigWig files. 

## Snakepipelines
Analysis pipelines built using Snakemake to output desired results (rna-seq counts, etc.)
Each folder contains at least a Snakefile (master file that tells how which output files you want), a configuration file (for custom-parameters that should not be in the main pipeline).
