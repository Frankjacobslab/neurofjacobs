##########################################################
## Number of detected repeats per class and per time point
##########################################################

#Load library
library(ggplot2)
library(RColorBrewer)
library(dplyr) # to perform split/apply/combine (like averaging)
library(stats)
library(optparse)
library(data.table)
library(reshape2) # especially for melt = wide to long format
library(gplots)
#library(multidplyr) # install.packages("devtools") and then devtools::install_github("hadley/multidplyr")

# melt (wide to long format) and add experimental timepoint per sample  
# load experimental design and add condition per sample
design = read.delim(file = "design.txt", header = T,stringsAsFactors = F)
scaled_counts_repeats = read.delim(file = "scaled.repeat.counts_>0reads_test.txt", header = T,stringsAsFactors = F)
head(scaled_counts_repeats)
m_counts = melt(scaled_counts_repeats,id.vars = c("repClass","ids"),value.name = "counts",variable.name = "sample",factorsAsStrings = F)
m_counts = left_join(x = m_counts,y = design[c("sample","condition")],by="sample")

# average sample counts that belong to same time point
m_counts_avg = m_counts %>%
  group_by(ids,condition,repClass) %>%
  summarise(counts = mean(counts))
write.table(m_counts_avg, file = "scaled_repeat_counts_>0reads_mean.txt",quote = F,row.names =F,sep="\t")

#filter out repeats that are expressed below a certain threshold
repeats_filtered = filter(m_counts_avg, D0>20 OR D28>20 OR D42>20)

# filter out repeats that are not expressed at any time point
#above_zeros = m_counts %>% 
  #group_by(ids) %>%
 # summarise(sum = sum(counts)) %>%
 # filter(sum > 0) 
#above_zeros = unique(above_zeros$ids)
#m_counts_avg_nonzeros = m_counts_avg[which(m_counts_avg$ids %in% above_zeros),]

# number of detection time for each repeat class
elementsDetectedPerClassAndTimepoint = m_counts_avg_nonzeros %>%
  group_by(repClass,time) %>%
  summarise(nb = sum(counts != 0))

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