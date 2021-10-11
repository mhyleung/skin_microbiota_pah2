#Generate bray-curtis and jaccard
library(readr)
library(dplyr)
library(tidyr)
library(wilkoxmisc)
library(vegan)
library(reshape2)
library(RVAideMemoire)
library(ggplot2)
#install.packages("BiocManager")
#BiocManager::install("phyloseq")

tax_cast <- read.table("../kraken2/bracken2_all_samples_wide_0_00001_adjusted_for_human_clean_no_con.txt",header=T,row.names=1)

tax_cast <- t(tax_cast)
distance_b <- vegdist(tax_cast, method="bray")
distance_j <- vegdist(tax_cast, method="jaccard")

BrayMatrix <- as.matrix(distance_b)

#Prediction strength
ps <- prediction.strength(BrayMatrix,M=100)
# print the prediction strength values
#Prediction strength 
#Clustering method:  kmeans 
#Maximum number of clusters:  10 
#Resampled data sets:  100 
#Mean pred.str. for numbers of clusters:  1 0.8076294 0.5627498 0.6749701 0.504831 0.448959 0.3756168 0.3318927 0.2797484 0.2374394 
#Cutoff value:  0.8 
#Largest number of clusters better than cutoff:  2

#Calculate number of clusters and their clustering strength (cluster package)
library(cluster)
cluster <- pam(BrayMatrix,2,cluster.only=TRUE)
write.tidy(cluster,"PAM_cutotype.txt")