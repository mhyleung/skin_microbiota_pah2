library(reshape2)
library(ggplot2)
library(ggsci)
library(wilkoxmisc)

#Open wide table and convert to long
Table <- read.tidy("bracken2_all_samples_long_percentage_wide_newer_top12.txt")

Table <- melt(Table, id.vars="Species")

#Add metadata
Meta <- read.tidy("../beta_adonis2/PCoA_coordinates_kraken2_w_ctypes.txt")[c("Sample","Type")]
Merge <- merge(Table,Meta,by.x="variable",by.y="Sample",all.x = TRUE)

Merge$Species <- factor(Merge$Species, levels = c("Cutibacterium acnes","Staphylococcus aureus","Micrococcus luteus","Aeromonas hydrophila","Moraxella osloensis","Burkholderia dolosa","Pseudomonas fluorescens","Cutibacterium granulosum","Kocuria palustris","Janibacter indicus","Escherichia coli","Kytococcus sedentarius","Minor/Unclassified"))

#Plot
Plot <- ggplot(Merge,aes(x=variable,y=value,fill=Species))
Plot <- Plot + geom_bar(stat="identity") + theme_classic() + scale_y_continuous(expand = c(0,0)) + facet_wrap(~Type,scales="free_x")
Plot <- Plot + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())
Plot <- Plot + scale_fill_igv()
ggsave("kraken2_bacteiral_top12_by_cutotype.pdf",width=8,height=6,units="in")