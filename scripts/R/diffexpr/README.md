# Rscripts directory
This directory contains R scripts (mostly Rmarkdown format) for various uses (data retrieval, differential expression, etc.)

## Differential expression

### deseq.Rmd
This script takes as input the __counts__ from RNA-Seq read alignments (output of snakepipelines/rnaseq Snakemake pipeline) and a __design__ file (experimental design) . It performs:
*  __Scaling__ (unproperly termed "normalization"): `estimateSizeFactors` function
*  __Variance estimation__: `estimateDispersions` function
* __Differential expression__: `nbinomWaldTest` 

### Input files format
The __counts__ file must have the following structure: 

 Geneid    | mDA_D0_2 | mDA_D0_3 | mDA_D0_1 | mDA_D28_3 | mDA_D28_2 | mDA_D28_1 | mDA_D42_1 | mDA_D42_2 | mDA_D42_3 
-----------|----------|----------|----------|-----------|-----------|-----------|-----------|-----------|-----------
 653635    | 97       | 84       | 71       | 1960      | 840       | 765       | 822       | 891       | 1157      
 100422834 | 0        | 0        | 0        | 0         | 0         | 0         | 0         | 0         | 0         
 645520    | 0        | 0        | 0        | 0         | 0         | 0         | 0         | 0         | 0         
 79501     | 0        | 0        | 0        | 0         | 0         | 0         | 0         | 0         | 0         
 729737    | 49       | 55       | 128      | 2973      | 1252      | 965       | 1910      | 1785      | 2514      
 100507658 | 0        | 0        | 0        | 1         | 0         | 0         | 0         | 0         | 2         
 100132287 | 2        | 8        | 12       | 557       | 188       | 138       | 170       | 168       | 240       

 Geneid    | mDA_D0_2 | mDA_D0_3 | mDA_D0_1 | mDA_D28_3 | mDA_D28_2 | mDA_D28_1 | mDA_D42_1 | mDA_D42_2 | mDA_D42_3 
-----------|----------|----------|----------|-----------|-----------|-----------|-----------|-----------|-----------
 653635    | 97       | 84       | 71       | 1960      | 840       | 765       | 822       | 891       | 1157      
 100422834 | 0        | 0        | 0        | 0         | 0         | 0         | 0         | 0         | 0         


The __design__ file must have the following structure:

 sample    | condition | libType    
-----------|-----------|------------
 mDA_D0_2  | D0        | paired-end 
 mDA_D0_3  | D0        | paired-end 
 mDA_D0_1  | D0        | paired-end 
 mDA_D28_3 | D28       | paired-end 
 mDA_D28_2 | D28       | paired-end 
 mDA_D28_1 | D28       | paired-end 
 mDA_D42_1 | D42       | paired-end 
 mDA_D42_2 | D42       | paired-end 
 mDA_D42_3 | D42       | paired-end 

### Running the script
1.  First move the first move R working directory to that of the `deseq.Rmd` script:
2.  Change the input files in the `params` header of the `deseq.Rmd` file. You need to specify a `counts` file and a `design` file. 
    + `counts`: `../../results/20160928
3.  Run `ezknitr::ezknit(file = "deseq.Rmd",out_dir = "../../analysis/20161003_diffexpression/",fig_dir = "figs")`	

