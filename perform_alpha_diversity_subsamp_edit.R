library(readr)
library(dplyr)
library(reshape2)
library(tidyr)

alpha <- read_tsv("merged_abundance_table_species_subsamp.txt") %>%
	gather(Variable, abundance, 2:(ncol(alpha))) %>%
	filter(abundance >0)
names(alpha)[1:2] <- c("Species", "Sample")

meta <- read_tsv("metadata_new.txt") %>%
	select(Sample, City)
	
merge <- left_join(alpha, meta)

richness <- merge %>%
	select(-abundance, -Species) %>%
	mutate(richness = 1) %>%
	group_by(Sample) %>%
	mutate(richness= sum(richness)) %>%
	ungroup() %>%
	unique()
write_tsv(richness,"richness_rarefied_global_samples_subsamp.txt")
	
# calculate the shannon diversity
library(vegan)
alpha[1:6,1:6]

alpha <- read.table("../bracken_subsamp/bracken_table_species_subsamp_transposed.txt", header=T, sep="\t", fill =TRUE, quote="", row.names = NULL)
Sample <- colnames(alpha[2:ncol(alpha)])
Species <- alpha[1:nrow(alpha),1]
names(alpha) <- NULL
alpha[,1] <- NULL
#transpose the data frame
alpha_t <- data.frame(t(alpha))
alpha_t[is.na(alpha_t)] <- 0
colnames(alpha_t) <- Species
#calculate the shannon diversity for global samples
Shannon <- diversity(alpha_t, "shannon")
#merge the sample info with alpha diversity
data <- data.frame(cbind(Sample,Shannon))
#load the meta table
meta <- read_tsv("../meta/meta_pahoh_new.txt")
#add meta info
merge <- left_join(data, meta)
write_tsv(merge, "shannon_diversity_rarefied_samples_bracken_subsamp.txt")

library(wilkoxmisc)
Table <- read.tidy("shannon_diversity_rarefied_samples_bracken_subsamp.txt")
meta <- read.tidy("../meta/meta_pahoh_cutotype.txt")
Table <- merge(Table,meta,by="Sample",all.X=TRUE)
Plot <- ggplot(Table, aes(x = City, y = Shannon, fill = City))
Plot <- Plot + geom_boxplot() + facet_wrap(~skin_imperfection_Para40)
Plot <- Plot + theme_classic()
Plot <- Plot + scale_fill_brewer(palette = "Dark2")
Plot <- Plot + theme(axis.title = element_text(size=14, face="bold"))
Plot <- Plot + ylab(paste0("Shannon Diversity")) + xlab(paste0("City"))
ggsave("shannon_city_acne_bracken.pdf",width=4,height=4,units="in")


#statistical test
#Kraken2
noacne <- Table[which(Table$skin_imperfection_Para40 == "No Acne") ,]
test <- wilcox.test(Shannon~City,data=noacne)
#data:  Shannon by City
#W = 545, p-value = 0.01554
#alternative hypothesis: true location shift is not equal to 0

acne <- Table[which(Table$skin_imperfection_Para40 == "Acne") ,]
test <- wilcox.test(Shannon~City,data=acne)
#data:  Shannon by City
#W = 679, p-value = 0.1303
#alternative hypothesis: true location shift is not equal to 0

test <- wilcox.test(Shannon.x~Cutotype,data=Table)
#Wilcoxon rank sum test with continuity correction
#data:  Shannon.x by Cutotype
#W = 75, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

noacne <- Table[which(Table$skin_imperfection_Para40.x == "No Acne") ,]
test <- wilcox.test(Shannon.x~Cutotype,data=noacne)
#W = 28, p-value = 3.021e-12
#alternative hypothesis: true location shift is not equal to 0

noacne <- Table[which(Table$skin_imperfection_Para40.x == "Acne") ,]
test <- wilcox.test(Shannon.x~Cutotype,data=noacne)
#W = 6, p-value = 9.88e-16
#alternative hypothesis: true location shift is not equal to 0
