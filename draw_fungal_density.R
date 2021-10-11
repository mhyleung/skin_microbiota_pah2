#Draw fungal density plot
library(wilkoxmisc)

Table <- read.tidy("findfungi_output_normalized.txt")[c("Sample", "Fungal_Abundance")]
Table_unique <- unique(Table)

Meta <- read.tidy("../meta/meta_sample_city_only.txt")

Merge <- merge(Table_unique,Meta, by="Sample",all.X=TRUE)

#Plot density plot
Plot <- ggplot(Merge,aes(x=Sample,y=Fungal_Abundance,fill=City))
Plot <- Plot + geom_bar(stat="identity") + theme_classic() + scale_y_continuous(expand = c(0,0)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
ggsave("fungal_overall_abundance.pdf",width=5,height=3,units="in")