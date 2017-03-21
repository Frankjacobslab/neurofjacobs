#!/usr/bin/env Rscript
#####################################################################
### Basic statistics and over-representation analysis for repeats 
# author: marc galland
# contact: m.galland@uva.nl
# version 1.0 created on December 14, 2016

# Three of the four inputs come from the Subread feature count analysis of RNA-Seq mapping files (BAM files)
# count_table should not be scaled/normalized in any way (will be by DESeq)

#Usage
# Rscript --vanilla repeats.R [path/to/count_table] [path/to/repeat_classification] [path/to/output_folder]

# to be added soon: [path/to/count_summary]

# Arguments: 
#   count_table: table of raw counts as given by Feature Counts program of Subread
#   count_summary: statistical summary of feature count analysis (number of assigned / non assigned)
#   repeat_classification: classification of repeats per family
#   output_folder: a path to a folder that will contain all results

#####################################################################

args = commandArgs(trailingOnly=TRUE)

setwd("~/SURFdrive/workspace_frank_jacobs_group/neurofjacobs/results/20161109_bigwigs_repeat_counts/")

# load library
library(ggplot2)
library(RColorBrewer)
library(dplyr) # to perform split/apply/combine (like averaging)
library(stats)
library(optparse)
library(data.table)
library(reshape2) # especially for melt = wide to long format
library(gplots)
library(multidplyr) # install.packages("devtools") and then devtools::install_github("hadley/multidplyr")
library(DESeq2)

##########################################
## load data and add information to counts
##########################################
count_table = "repeat_counts.txt"
#count_summary = "repeat_counts.txt.summary"
repeat_classification = "../../data/02.refs/repeat_classification.txt"
xp_design = "../../data/01.seqruns/20160912_dopamine_spheres/design.txt"
outdir = "repeat_analysis/" 

# load repeat counts
# remove unnecessary columns (Chr Start End etc.)
counts = fread(count_table,header = T,stringsAsFactors = F)

# load mapping summary 
#count_summary = read.delim(count_summary,header = T,stringsAsFactors = F)

# load repeat classification (units to class)
classif = read.delim(repeat_classification,header=T,stringsAsFactors = F)

##########################################################
## Add repeat name and create new repeat id
##########################################################

# rename first column to match that of repeat classification
# merge to get the type of each repeat
colnames(counts)[1] = "repName"
counts = left_join(x = counts,y = classif,by="repName")

# create a unique id (because repeats have....repeated names!)
counts$ids = with(data = counts,paste(repName,Chr,Start,End,Strand,sep = "_"))
counts = subset(counts,select = -c(1,2,3,4,5,6))


##############
# Scaled counts
##############
size_factors = read.delim("sizeFactors_repeats.txt",header=T,stringsAsFactors = F)
design = read.delim(xp_design,header=T,stringsAsFactors = F)

# create count matrix
deseq.matrix = counts
deseq.matrix$ids <- NULL
deseq.matrix$repClass <- NULL

# convert to DESeq count dataset
bob = DESeqDataSetFromMatrix(countData=as.matrix(deseq.matrix),colData=design,design=~time)
bob = estimateSizeFactors(bob)
scaled.repeat.counts=counts(bob,normalized=T)

# add rep name and class
scaled.repeat.counts=as.data.frame(scaled.repeat.counts)
scaled.repeat.counts$ids = counts$ids
scaled.repeat.counts$repClass = counts$repClass

# add a column that contains sum of row
# filter lines where sum of scaled counts == 0
sumOfRows = rowSums(as.matrix(scaled.repeat.counts[,1:9]))
scaled.repeat.counts = scaled.repeat.counts[which(sumOfRows > 0),]

# add back the repeat name as a column
colsplitted = colsplit(scaled.repeat.counts$ids,pattern = "_",names = c("repName","chr","start","end","strand"))
final = cbind.data.frame(colsplitted,scaled.repeat.counts) 

# write to table
write.table(final,file = file.path(outdir,"scaled.repeat.counts.txt"),quote = F,row.names =F,sep="\t") 


##########################################################
## Number of detected repeats per class and per time point
##########################################################

# melt (wide to long format) and add experimental timepoint per sample  
# load experimental design and add condition per sample
m_counts = melt(scaled.repeat.counts,id.vars = c("repClass","ids"),value.name = "counts",variable.name = "sample",factorsAsStrings = F)
m_counts = left_join(x = m_counts,y = design[c("sample","time")],by="sample")

# average sample counts that belong to same time point
m_counts_avg = m_counts %>%
    group_by(ids,time,repClass) %>%
    summarise(counts = mean(counts))

# filter out repeats that are not expressed at any time point
above_zeros = m_counts %>% 
    group_by(ids) %>%
    summarise(sum = sum(counts)) %>%
    filter(sum > 0) 
above_zeros = unique(above_zeros$ids)
m_counts_avg_nonzeros = m_counts_avg[which(m_counts_avg$ids %in% above_zeros),]

# number of detection time for each repeat class
elementsDetectedPerClassAndTimepoint = m_counts_avg_nonzeros %>%
  group_by(repClass,time) %>%
  summarise(nb = sum(counts != 0))


###################################################
## Percentages of reads mapped to each repeat class
###################################################

# import mapping logs and merge them altogether
mappings = list.files(path="mappinglogs",full.names = T)
mappings = lapply(mappings,FUN = function(x) {read.delim(x,header = F,sep = "\t",stringsAsFactors = F)})
names(mappings) = design$sample
for (i in seq_along(mappings)){
  newColname = names(mappings)[i]
  colnames(mappings[[i]]) = c("params",newColname)
}
all.mappings = Reduce(function(...) merge(...,by="params",sort=F),mappings)

# melt data frame
m_all.mappings = melt(data = all.mappings,id.vars = c("params"),factorsAsStrings=T,variable.name = "sample")

# number of trimmed reads used as input for STAR mapping 
# these numbers will be used to compute percentages of reads mapped to each repeat category
inputReads = m_all.mappings[grep(pattern = "Number of input reads",x = m_all.mappings$params),]
colnames(inputReads)[3] = "input_reads"
inputReads = inputReads[order(inputReads$sample),]

# add number of input reads to count table
# sum per sample and per repeat class
sumClassCounts = m_counts %>%
  group_by(repClass,sample) %>% 
  summarise(sum = sum(counts))

sumClassCounts = left_join(sumClassCounts,inputReads,by="sample")
sumClassCounts$percents = sumClassCounts$sum / as.numeric(sumClassCounts$input_reads) * 100
write.table(sumClassCounts,file = file.path(outdir,"percentages_mapped_per_class.txt"),sep = "\t",row.names = F,quote = F)


########
## Plots
########

# color blind compatible
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# define a custom theme
my_theme <- theme_grey() + 
  theme(
    axis.text = element_text(face = "bold",colour = "black",size = 8)
  )

# plot quantiles from counts
quantiles = data.frame(q = quantile(m_counts_avg_nonzeros$counts,c(seq(from=0.1,to=0.9,by=0.1),seq(from=0.91,to=0.99,by=0.01))))
quantiles$x.axis = c(
  paste(seq(10,90,10),"%",sep=""),
  paste(seq(91,99,1),"%",sep=""))

g <- ggplot(data = quantiles,aes(x = x.axis,y=q)) +
  geom_bar(stat="identity") +
  labs(x = "Quantiles",y="Scaled counts")
print(g)
ggsave(filename = file.path(outdir,"quantiles.png"),plot = g,width = 7,height = 5,dpi = 400)

# Plots of the number of repeats detected (scaled counts > 0) for each family and at each time point
p <- ggplot(data = elementsDetectedPerClassAndTimepoint,mapping = aes(x = time,y = nb,fill=time)) +
  geom_bar(stat="identity") +
  facet_wrap(~ repClass,nrow=4,scales="free") +
  labs(x = "",y = "Number of repeat elements detected") +
  my_theme +
  theme(
    axis.text.x=element_blank(), # removes x axis
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 6)
  ) +
  scale_fill_manual(values = cbPalette)
print(p)
ggsave(filename = file.path(outdir,"n_repeats_per_timepoint.png"),plot = p,width = 7,height = 5,dpi = 600)
ggsave(filename = file.path(outdir,"n_repeats_per_timepoint.svg"),plot = p,width = 7,height = 5)

# plot percentages of mapped reads per repeat class
# add time for plots
sumClassCounts = left_join(sumClassCounts,design,by="sample")
p2 <- ggplot(data = sumClassCounts,mapping = aes(x = sample,y = percents,fill=time)) +
  geom_bar(stat="identity") +
  facet_wrap(~ repClass,nrow=4,scales="free") +
  labs(x = "",y = "Fraction of trimmed reads that map to repeat family (%)") +
  my_theme +
  theme(
    axis.text.x=element_blank(), # removes x axis
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 6)
  ) +
  scale_fill_manual(values=cbPalette)
print(p2)
ggsave(filename = file.path(outdir,"percent_mapped_per_repeat_class.png"),plot = p2,width = 7,height = 5,dpi = 600)
ggsave(filename = file.path(outdir,"percent_mapped_per_repeat_class.svg"),plot = p2,width = 7,height = 5)

sessionInfo()
