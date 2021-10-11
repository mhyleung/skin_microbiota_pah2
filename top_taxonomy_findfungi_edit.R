library(wilkoxmisc)
#Open findfungi file
Table <- read.tidy("findfungi_output_normalized.txt")[c("Sample","Species","RelativeAbundnce_communuty")]

Table <- melt(Table)

Species <- ddply(Table, .(Sample, Species), summarise, value = sum(value))
names(Species)[3] <- "RelativeAbundance"

Species <- collapse.taxon.table(Species, n = 12, Rank = "Species")

Table2 <- read.tidy("findfungi_output.txt")[c("Sample","Overall")]
Table2 <- unique(Table2)
Species <- merge(Species,Table2,by="Sample",all.X=TRUE)
write.tidy(Species,"top12findfungi.txt")


Species <- read.tidy("top12findfungi.txt")
Species$Species <- factor(Species$Species,levels=c("Malassezia restricta","Malassezia globosa","Alternaria alternata","Preussia sp. BSL10","Puccinia arachidis","Alternaria arborescens","Chrysosporium queenslandicum",
                                                    "Amauroascus mutatus","Amauroascus niger","Clavaria fumosa","Malassezia furfur","Minor/Unclassified"))

Plot <- ggplot(Species,aes(x=Sample))
Plot <- Plot + geom_bar(aes(y=RelativeAbundance,fill=Species),stat="identity") + facet_grid(~City,scales="free") + scale_y_continuous(expand = c(0,0)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_fill_brewer(palette="Paired")
ggsave("findfungi_top12.pdf",width=8,height=6,units="in")
