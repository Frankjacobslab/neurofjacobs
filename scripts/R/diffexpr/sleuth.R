#!/usr/bin/env Rscript
#####################################################################
### Sleuth analysis to call differential expression on repeats
# author: marc galland
# contact: m.galland@uva.nl
# version 1.0 created on January 27, 2017

# Three of the four inputs come from the Subread feature count analysis of RNA-Seq mapping files (BAM files)
# count_table should not be scaled/normalized in any way (will be by DESeq)

#Usage
# Rscript --vanilla repeats.R [datadir] [info] [outdir]
# Arguments: 
#     datadir: a path to a directory containing Kallisto results (one file per sample)
#     info: experimental design (condition to samples)
#     outdir: a path to a folder that will contain all results
#####################################################################

## load librairies
## allow multi-threading
library(sleuth)
library(dplyr)
options(mc.cores = 4L)

##################
## Parse arguments
##################
args = commandArgs(trailingOnly = TRUE)
datadir = args[1]
infoFile = args[2]
outdir = args[3]

#####################
## Load kallisto data
#####################
# Configurate file for Sleuth
sample_id = dir(datadir)
kal_dirs <- sapply(sample_id, function(id) {file.path(datadir, id)})
s2c = read.table(infoFile,header=T,stringsAsFactors=FALSE)
s2c = dplyr::mutate(s2c,path=kal_dirs)

# Creation of the sleuth object 
so = sleuth_prep(s2c, ~ time)

# Extraction of the normalized expression data 
df = kallisto_table(so,normalized = T,use_filtered = T)

# Sleuth analysis
so = sleuth_fit(so) # full model
so = sleuth_fit(so,formula = ~ 1,fit_name = "reduced")
so = sleuth_lrt(so,"reduced","full")

#sleuth_results(likelihood ratio test)
sr <- sleuth_results(so,"reduced:full","lrt")
write.table(sr,file.path(outdir,"sleuth_results.txt",sep="\t",row.names=F)

