library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(wilkoxmisc)

#identify core pathways
table <- read.tidy("joined_pathwayabundance_relab_stratified_fixed.tsv")

table <- table %>%
gather(Sample, Abundance, 3:(ncol(table))) %>%
    dplyr::filter(Abundance > 0)

table <- table %>%
    mutate(Sample = gsub( "_"," ",Sample)) %>%
    mutate(Sample = word(Sample,1))

meta <- read_tsv("meta_sample_city_cutotype.txt") %>%
    select(Sample, Type)

table <- left_join(table, meta)

table <- table %>%
    mutate(Species = ifelse(Species == "unclassified","unclassified", "classified"))

table <- table %>% dplyr::filter(Type == "C2")

table <- table %>%
    group_by(Pathway, Species, Sample) %>%
    mutate(Abundance=sum(Abundance)) %>%
    ungroup() %>%
    unique()

table <- table %>%
    group_by(Pathway, Sample) %>%
    mutate(Total_abundance = sum(Abundance)) %>%
    ungroup() %>%
    unique()

table <- table %>%
    mutate(Proportion = (Abundance/Total_abundance)*100)
length(unique(table$Pathway))
#245 for C1, 298 for C2
length(unique(table$Sample))
#73 for C1, 49 for C2

# remove pathways with mean contribution from unclassified groups greater than 25%
unclassified <- table %>% dplyr::filter(Species == "unclassified" & Proportion > 25)
length(unique(unclassified$Pathway))
#0 for C1 and C2

table_filtered <- table %>%
    dplyr::filter(!Pathway %in% unclassified$Pathway)

table_filtered2 <- table_filtered %>%
    dplyr::filter(Species == "classified") %>%
    group_by(Pathway) %>%
    dplyr::summarise(number= n())

merge <- left_join(table_filtered, table_filtered2) %>%
    mutate(Prevalence = (number/49)*100) %>% #change to total number of samples in particular group you are looking at
    dplyr::filter(Prevalence > 75) %>%
    select(Pathway, Prevalence) %>%
    unique()

write_tsv(merge, "core_pathways_C2_new.txt")

library(readr)
library(dplyr)
library(tidyr)
#renormalize pathway abundance table by excluding unclassified groups
table <- read_tsv("humann_pathabundance_relab_stratified_C2.txt") %>%
    dplyr::filter(! Species == "unclassified")

paw <- table %>% dplyr::filter(Pathway == "VALSYN-PWY: L-valine biosynthesis")
#paw <- paw[,-1]

paw <- paw %>%
gather(Sample, Abundance, 3:(ncol(paw))) %>%
dplyr::filter(Abundance > 0) %>%
select(Sample, Species, Abundance)

#calculate Gini-Simpson index
paw <- as.data.frame(paw)
library(diverse)
diversity <- diverse::diversity(paw, type = "gini-simpson")
mean(diversity$gini.simpson)

library(wilkoxmisc)
library(reshape2)
paw_cast <- dcast(paw, Sample ~ Species, value.var = "Abundance", fill = 0)
write_tsv(paw_cast, "paw_cast_C1.txt")


paw_cast <- read.table("paw_cast_C1.txt",header=T,row.names=1)
#distance <- distance(paw_cast, method = "bray-curtis") in R package "ecodist
library(vegan)
distance <- vegdist(paw_cast, method="bray")
mean(distance)

