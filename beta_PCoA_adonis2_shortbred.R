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

tax_cast <- read.table("combined_shortbred_output_wide_transposed.txt",header=T,row.names=1)

tax_cast <- t(tax_cast)
distance_b <- vegdist(tax_cast, method="bray")
distance_j <- vegdist(tax_cast, method="jaccard")

BrayMatrix <- as.matrix(distance_j)


write.table(BrayMatrix,"BC_resistome_matrix.txt",row.names=TRUE)
Meta <- read_tsv("../meta/meta_pahoh_cutotype.txt")[c("Sample","City","Cutotype_PAM","Shannon","Age_Group","Cutibacterium_granulosum","Saureus","Micrococcus_luteus","Aeromonas_hydrophila","Moraxella_osloensis","Pseudomonas_aeruginosa","skin_imperfection_Para40","Face_Pigmentation_Frequency","Nicotine","Cotinine","Acenaphtylene","Acenaphtene","Fluorene",
                                                 "Phenanthrene","Anthracene","Fluoranthene","Pyrene","Benzo[a]anthracene","Chrysene","Benzo[b]fluoranthene","Benzo[k]fluoranthene","Benzo[a]pyrene","Indeno[1,2,3-cd]pyrene","Dibenzo[a,h]anthracene","Benzo[ghi]perylene")]

PCoA <- cmdscale(BrayMatrix, k = 2, eig = TRUE)
DF <- data.frame(Sample = row.names(PCoA$points), PCoA1 = PCoA$points[,1], PCoA2 = PCoA$points[,2], row.names = NULL)
Eigenvalues <- eigenvals(PCoA)
VarianceExplained <- Eigenvalues /sum(Eigenvalues)
VarianceExplained1 <- 100 * signif(VarianceExplained[1], 2)
VarianceExplained2 <- 100 * signif(VarianceExplained[2], 2)
PCoA <- merge(DF, Meta, by = "Sample", all.x = TRUE)
write.tidy(PCoA,"cutotype_PAM_PCoA_coords.txt")

AllSamples <-data.frame(Sample = row.names(BrayMatrix))
AllSamples <- merge(AllSamples, Meta, by = "Sample", all.x = TRUE)
#Check that UniFrac matrix rows match samples table (for ANOSIM)
sum(row.names(BrayMatrix)==AllSamples$Sample) == length(AllSamples$Sample)

BrayCurtis <- as.dist(distance_b)
Jaccard <- as.dist(distance_j)

xlab <- paste0("PCoA1 (Variance Explained: 29%)") #Check VarianceExplained1
ylab <- paste0("PCoA2 (Variance Explained: 15%)") #Check VarianceExplained2

#PCoA <- read.tidy("cutotype_PAM_PCoA_coords.txt")
Plot <- ggplot(PCoA, aes(x = PCoA1, y = PCoA2))
Plot <- Plot + geom_point(aes(colour=Cutotype_PAM)) + stat_ellipse(geom = "polygon",alpha = 0.09, aes(fill = Cutotype_PAM))
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
Plot <- Plot + theme_classic()
#Plot <- Plot + scale_color_brewer(palette = "Set2")
Plot <- Plot + theme(axis.title=element_text(size=14)) + scale_size_continuous(range=c(0.1,5))
Plot <- Plot + theme(axis.text=element_text(size=12))
Plot <- Plot + theme(legend.title = element_blank(), legend.position="right")
Plot <- Plot + theme(legend.text = element_text(size=14, face = "bold"))
ggsave("shortbred_bray_curtis_by_cutotype.pdf")
write.table(PCoA,"PCoA_coordinates_kraken2.txt")

#Perform Adonis(permanova)
adonis(formula = BrayCurtis ~ Cutotype_PAM*skin_imperfection_Para40*City*Age_Group*Face_Pigmentation_Frequency, data = AllSamples, permutations = 999)
Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Cutotype_PAM                                                             1    4.5078  4.5078 25.5098 0.17713  0.001 ***
skin_imperfection_Para40                                                 1    0.1341  0.1341  0.7589 0.00527  0.644    
City                                                                     1    0.1506  0.1506  0.8523 0.00592  0.528    
Age_Group                                                                3    0.4036  0.1345  0.7613 0.01586  0.779    
Face_Pigmentation_Frequency                                              1    0.1591  0.1591  0.9005 0.00625  0.510    
Cutotype_PAM:skin_imperfection_Para40                                    1    0.1502  0.1502  0.8502 0.00590  0.566    
Cutotype_PAM:City                                                        1    0.1276  0.1276  0.7224 0.00502  0.677    
skin_imperfection_Para40:City                                            1    0.1111  0.1111  0.6287 0.00437  0.773    
Cutotype_PAM:Age_Group                                                   3    0.4272  0.1424  0.8059 0.01679  0.732    
skin_imperfection_Para40:Age_Group                                       3    0.4799  0.1600  0.9052 0.01886  0.583    
City:Age_Group                                                           3    0.6384  0.2128  1.2043 0.02509  0.217    
Cutotype_PAM:Face_Pigmentation_Frequency                                 1    0.2466  0.2466  1.3957 0.00969  0.157    
skin_imperfection_Para40:Face_Pigmentation_Frequency                     1    0.1153  0.1153  0.6523 0.00453  0.757    
City:Face_Pigmentation_Frequency                                         1    0.1772  0.1772  1.0025 0.00696  0.403    
Age_Group:Face_Pigmentation_Frequency                                    3    0.5420  0.1807  1.0224 0.02130  0.437    
Cutotype_PAM:skin_imperfection_Para40:City                               1    0.1467  0.1467  0.8304 0.00577  0.572    
Cutotype_PAM:skin_imperfection_Para40:Age_Group                          3    0.4005  0.1335  0.7555 0.01574  0.785    
Cutotype_PAM:City:Age_Group                                              3    0.6263  0.2088  1.1813 0.02461  0.273    
skin_imperfection_Para40:City:Age_Group                                  3    0.4338  0.1446  0.8183 0.01705  0.703    
Cutotype_PAM:skin_imperfection_Para40:Face_Pigmentation_Frequency        1    0.1244  0.1244  0.7040 0.00489  0.678    
Cutotype_PAM:City:Face_Pigmentation_Frequency                            1    0.1635  0.1635  0.9254 0.00643  0.470    
skin_imperfection_Para40:City:Face_Pigmentation_Frequency                1    0.1085  0.1085  0.6143 0.00427  0.815    
Cutotype_PAM:Age_Group:Face_Pigmentation_Frequency                       3    0.4440  0.1480  0.8375 0.01745  0.688    
skin_imperfection_Para40:Age_Group:Face_Pigmentation_Frequency           1    0.1393  0.1393  0.7882 0.00547  0.596    
City:Age_Group:Face_Pigmentation_Frequency                               2    0.3900  0.1950  1.1035 0.01532  0.350    
Cutotype_PAM:skin_imperfection_Para40:City:Age_Group                     2    0.4823  0.2411  1.3645 0.01895  0.147    
Cutotype_PAM:skin_imperfection_Para40:City:Face_Pigmentation_Frequency   1    0.1863  0.1863  1.0542 0.00732  0.390    
Cutotype_PAM:City:Age_Group:Face_Pigmentation_Frequency                  1    0.1800  0.1800  1.0188 0.00707  0.411    
Residuals                                                               75   13.2533  0.1767         0.52076           
Total                                                                  123   25.4498                 1.00000       

#Classify markers into classes
Table <- read.tidy("combined_shortbred_output_wide_family_level.txt")
Table <- melt(Table)
aro <- read.tidy("aro.tsv")
Merge <- merge(aro,Table,by.x="Accession",by.y="Family", all.X=TRUE)

#Merge with cutotype data
Merge <- merge(Merge,Meta,by.x="variable",by.y="Sample",all.X+TRUE)
write.tidy(Merge,"ARG_family_w_cutotype.txt")

#Generate plot for ShortBRED
#Draw plots for significant maaslin correlations for PAH
Table <- read.tidy("ARG_family_w_cutotype.txt")
ARO_3000478 <- Table[which(Table$Accession == "ARO_3000478") ,]
Plot1 <- ggplot(ARO_3000478,aes(x=value,y=Dibenzo_a_h_anthracene,colour=Cutotype_PAM))
Plot1 <- Plot1 + geom_point()
