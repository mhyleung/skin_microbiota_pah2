library(readr)
library(dplyr)
library(tidyr)
library(readr)
library(dplyr)
library(tidyr)
library(wilkoxmisc)

table <- read_tsv("humann_pathabundance_relab_stratified_C2.txt") %>%
    dplyr::filter(!Species == "unclassified")
   
#use core pathway file from "identify_core_pathway" script to calculate gini and bray one by one
path <- table %>% dplyr::filter(Pathway == "VALSYN-PWY: L-valine biosynthesis")
path <- path[,-1]

path <- path %>%
gather(Sample, Abundance, 3:(ncol(path))) %>%
dplyr::filter(Abundance > 0) %>%
select(Sample, Species, Abundance)

#calculate Gini-Simpson index
path <- as.data.frame(path)
library(diverse)
diversity <- diverse::diversity(path, type = "gini-simpson")
mean(diversity$gini.simpson)

library(wilkoxmisc)
library(reshape2)
path_cast <- dcast(path, Sample ~ Species, value.var = "Abundance", fill = 0)
write_tsv(path_cast, "path_cast.txt")


path_cast <- read.table("path_cast.txt",header=T,row.names=1)
#distance <- distance(path_cast, method = "bray-curtis") in R package "ecodist
library(vegan)
distance <- vegdist(path_cast, method="bray")
mean(distance)

# make the plot of contributional diversity
library(readr)
library(dplyr)
library(tidyr)
library(readr)
library(dplyr)
library(tidyr)
table <- read_tsv("contributional_diversity_core_pathway_per_cutotype_new.txt")
table <- table %>%
    mutate(within_diversity = ifelse(gini.simpson < 0.5, "low", "high"))%>%
    mutate(between_diversity= ifelse(bray_curtis < 0.5, "low", "high"))

table2 <- table %>%
    mutate(group = ifelse(within_diversity == "low" & between_diversity == "low", "simple_conserved", NA)) %>%
    mutate(group = ifelse(within_diversity == "high" & between_diversity == "high", "complex_variable", group)) %>%
    mutate(group = ifelse(within_diversity == "high" & between_diversity == "low", "complex_conserved", group)) %>%
    mutate(group = ifelse(within_diversity == "low" & between_diversity == "high", "simple_variable", group))
write_tsv(table2,"contributional_diversity_core_pathway_per_cutotype_new.txt")

library(ggplot2)
table2 <- read_tsv("contributional_diversity_core_pathway_per_cutotype_new.txt")
Plot <- ggplot(table2, aes(x=gini.simpson, y=bray_curtis, color = group,shape=Shared))
Plot <- Plot + geom_point(size=2)
Plot <- Plot + theme_bw()
Plot <- Plot + xlab(paste0("Within-sample contributional diversity \n(Gini-Simpson)")) + ylab(paste0("Between-sample contributional diversity \n(Bray-Curtis)"))
Plot <- Plot + scale_color_brewer(palette = "Dark2")
Plot <- Plot + facet_wrap(~Cutotype)
Plot <- Plot + expand_limits(x = 0, y = 0)
Plot <- Plot + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),limits = c(0,1))
Plot <- Plot + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),limits = c(0,1))
Plot <- Plot + theme(axis.title = element_text(size = 14, face = "bold"))
Plot <- Plot + theme(legend.text=element_text(size = 12))
Plot <- Plot + theme(legend.title=element_text(size=14, face="bold")) + theme(legend.title = element_blank())
Plot <- Plot + geom_hline(yintercept = 0.5, colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed")
ggsave("contributional_diversity_new.pdf",width=6, height=5,units="in")

