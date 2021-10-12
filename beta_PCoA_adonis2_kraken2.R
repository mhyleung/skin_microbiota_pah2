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

write.table(BrayMatrix,"BC_tax_matrix.txt",row.names=TRUE)
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

xlab <- paste0("PCoA1 (Variance Explained: 44%)") #Check VarianceExplained1
ylab <- paste0("PCoA2 (Variance Explained: 9.1%)") #Check VarianceExplained2

PCoA <- read.tidy("cutotype_PAM_PCoA_coords.txt")
Plot <- ggplot(PCoA, aes(x = PCoA1, y = PCoA2))
Plot <- Plot + geom_point(aes(colour=Cutotype_PAM,size=Shannon)) + stat_ellipse(geom = "polygon",alpha = 0.09, aes(fill = Cutotype_PAM))
Plot <- Plot + xlab(xlab)
Plot <- Plot + ylab(ylab)
Plot <- Plot + theme_classic()
#Plot <- Plot + scale_color_brewer(palette = "Set2")
Plot <- Plot + theme(axis.title=element_text(size=14)) + scale_size_continuous(range=c(0.1,5))
Plot <- Plot + theme(axis.text=element_text(size=12))
Plot <- Plot + theme(legend.title = element_blank(), legend.position="right")
Plot <- Plot + theme(legend.text = element_text(size=14, face = "bold"))
ggsave("kraken2_bray_curtis_by_city_by_cutotype.pdf")
write.table(PCoA,"PCoA_coordinates_kraken2.txt")

#Perform Adonis(permanova)
adonis(formula = BrayCurtis ~ Cutotype_PAM*skin_imperfection_Para40*City*Age_Group*Face_Pigmentation_Frequency, data = AllSamples, permutations = 999)
                                                                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
skin_imperfection_Para40                                              1    0.1693  0.1693   2.020 0.01026  0.068 .  
City                                                                  1    0.4914  0.4914   5.862 0.02978  0.001 ***
Type                                                                  1    5.3306  5.3306  63.600 0.32303  0.001 ***
Age_Group                                                             3    0.1857  0.0619   0.739 0.01125  0.723    
Face_Pigmentation_Frequency                                           1    0.0478  0.0478   0.570 0.00289  0.756    
skin_imperfection_Para40:City                                         1    0.1355  0.1355   1.617 0.00821  0.141    
skin_imperfection_Para40:Type                                         1    0.0646  0.0646   0.770 0.00391  0.537    
City:Type                                                             1    0.1632  0.1632   1.947 0.00989  0.095 .  
skin_imperfection_Para40:Age_Group                                    3    0.2872  0.0957   1.142 0.01740  0.299    
City:Age_Group                                                        3    0.2170  0.0723   0.863 0.01315  0.565    
Type:Age_Group                                                        3    0.2662  0.0887   1.059 0.01613  0.353    
skin_imperfection_Para40:Face_Pigmentation_Frequency                  1    0.0377  0.0377   0.449 0.00228  0.883    
City:Face_Pigmentation_Frequency                                      1    0.1406  0.1406   1.677 0.00852  0.119    
Type:Face_Pigmentation_Frequency                                      1    0.0452  0.0452   0.539 0.00274  0.787    
Age_Group:Face_Pigmentation_Frequency                                 3    0.2206  0.0735   0.877 0.01337  0.594    
skin_imperfection_Para40:City:Type                                    1    0.1936  0.1936   2.310 0.01173  0.045 *  
skin_imperfection_Para40:City:Age_Group                               3    0.2647  0.0882   1.053 0.01604  0.405    
skin_imperfection_Para40:Type:Age_Group                               3    0.3025  0.1008   1.203 0.01833  0.250    
City:Type:Age_Group                                                   3    0.3984  0.1328   1.585 0.02414  0.102    
skin_imperfection_Para40:City:Face_Pigmentation_Frequency             1    0.1393  0.1393   1.662 0.00844  0.120    
skin_imperfection_Para40:Type:Face_Pigmentation_Frequency             1    0.1337  0.1337   1.596 0.00810  0.161    
City:Type:Face_Pigmentation_Frequency                                 1    0.0872  0.0872   1.040 0.00528  0.346    
skin_imperfection_Para40:Age_Group:Face_Pigmentation_Frequency        2    0.1037  0.0518   0.619 0.00628  0.798    
City:Age_Group:Face_Pigmentation_Frequency                            2    0.2651  0.1326   1.582 0.01607  0.112    
Type:Age_Group:Face_Pigmentation_Frequency                            2    0.1518  0.0759   0.905 0.00920  0.485    
skin_imperfection_Para40:City:Type:Age_Group                          1    0.1029  0.1029   1.227 0.00623  0.236    
skin_imperfection_Para40:City:Type:Face_Pigmentation_Frequency        1    0.1461  0.1461   1.743 0.00885  0.122    
skin_imperfection_Para40:City:Age_Group:Face_Pigmentation_Frequency   1    0.0406  0.0406   0.484 0.00246  0.826    
Residuals                                                            76    6.3700  0.0838         0.38601           
Total                                                               123   16.5021                 1.00000        

#Correlate PC coordinates and PAH levels
PCoA <- read.tidy("PCoA_coordinates_kraken2.txt")
#Start correlating PCoA coordinates and PAH exposure
cor.test(PCoA$PCoA1,PCoA$Fluoranthene,method="spearman")
Spearmans rank correlation rho

data:  PCoA$PCoA1 and PCoA$Fluoranthene
S = 241473, p-value = 0.01387
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.2213677

cor.test(PCoA$PCoA1,PCoA$Pyrene,method="spearman")
Spearmans rank correlation rho

data:  PCoA$PCoA1 and PCoA$Pyrene
S = 224716, p-value = 0.002119
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.2753995 

cor.test(PCoA$PCoA1,PCoA$Benzoaanthracene,method="spearman")
data:  PCoA$PCoA1 and PCoA$Benzoaanthracene
S = 232537, p-value = 0.005257
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.2501814 

cor.test(PCoA$PCoA1,PCoA$Chrysene,method="spearman")
data:  PCoA$PCoA1 and PCoA$Chrysene
S = 213292, p-value = 0.0004384
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.3122352

cor.test(PCoA$PCoA1,PCoA$Benzobfluoranthene,method="spearman")
data:  PCoA$PCoA1 and PCoA$Benzobfluoranthene
S = 243673, p-value = 0.01732
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.2142729

cor.test(PCoA$PCoA1,PCoA$Benzoghiperylene,method="spearman")
data:  PCoA$PCoA1 and PCoA$Benzoghiperylene
S = 232328, p-value = 0.005132
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.2508558 

cor.test(PCoA$PCoA2,PCoA$Acenaphtene,method="spearman")
data:  PCoA$PCoA2 and PCoA$Acenaphtene
S = 378060, p-value = 0.01492
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
-0.2190614 


#Draw correlation plots
PCoA <- read.tidy("PCoA_coordinates_kraken2.txt")

Plot <- ggplot(PCoA,aes(x=PCoA1,y=Fluoranthene,color=Cutotype))
Plot <- Plot + geom_point(size=2) + theme_classic()
Plot <- Plot + ylab(paste0("Fluoranthene (pg/mg hair)")) + xlab(paste0("PCoA1")) + theme(legend.position="none")
ggsave("PCoA1_fluoranthene_cutotype.pdf",width=3,height=3,units="in")

Plot <- ggplot(PCoA,aes(x=PCoA1,y=Pyrene,color=Cutotype))
Plot <- Plot + geom_point(size=2) + theme_classic()
Plot <- Plot + ylab(paste0("Pyrene (pg/mg hair)")) + xlab(paste0("PCoA1")) + theme(legend.position="none")
ggsave("PCoA1_pyrene_cutotype.pdf",width=3,height=3,units="in")

Plot <- ggplot(PCoA,aes(x=PCoA1,y=Benzo_a_anthracene,color=Cutotype))
Plot <- Plot + geom_point(size=2) + theme_classic()
Plot <- Plot + ylab(paste0("Benzo[a]anthracene (pg/mg hair)")) + xlab(paste0("PCoA1")) + theme(legend.position="none")
ggsave("PCoA1_benzoaanthracene_cutotype.pdf",width=3,height=3,units="in")

Plot <- ggplot(PCoA,aes(x=PCoA1,y=Chrysene,color=Cutotype))
Plot <- Plot + geom_point(size=2) + theme_classic()
Plot <- Plot + ylab(paste0("Chrysene (pg/mg hair)")) + xlab(paste0("PCoA1")) + theme(legend.position="none")
ggsave("PCoA1_chrysene_cutotype.pdf",width=3,height=3,units="in")

Plot <- ggplot(PCoA,aes(x=PCoA1,y=Benzo_b_fluoranthene,color=Cutotype))
Plot <- Plot + geom_point(size=2) + theme_classic()
Plot <- Plot + ylab(paste0("Benzo[b]fluoranthene (pg/mg hair)")) + xlab(paste0("PCoA1")) + theme(legend.position="none")
ggsave("PCoA1_Benzobfluoranthene_cutotype.pdf",width=3,height=3,units="in")

Plot <- ggplot(PCoA,aes(x=PCoA1,y=Benzo_ghi_perylene,color=Cutotype))
Plot <- Plot + geom_point(size=2) + theme_classic()
Plot <- Plot + ylab(paste0("Benzo[ghi]perylene (pg/mg hair)")) + xlab(paste0("PCoA1")) + theme(legend.position="none")
ggsave("PCoA1_Benzoghiperylene_cutotype.pdf",width=3,height=3,units="in")


#City-specific correlation analysis
PCoA <- read.tidy("PCoA_coordinates_kraken2.txt")
City <- PCoA[which(PCoA$City == "Baoding") ,]
cor.test(City$PCoA1,City$Pyrene,method="spearman")
data:  City$PCoA1 and City$Pyrene
S = 23664, p-value = 0.0178
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.3084746 

City <- PCoA[which(PCoA$City == "Dalian") ,]
cor.test(City$PCoA1,City$Benzoghiperylene,method="spearman")
Spearmans rank correlation rho

data:  City$PCoA1 and City$Benzoghiperylene
S = 32697, p-value = 0.04505
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.2514404 

#Perform dbRDA
library(wilkoxmisc)
BrayCurtis <- read.dist("BC_tax_matrix.txt")
Meta <- read.tidy("../meta/meta_pahoh_new_nospecialcharacter.txt")[c("Sample","City","Shannon","Age_Group","Cutibacterium_granulosum","Saureus","Micrococcus_luteus","Aeromonas_hydrophila","Moraxella_osloensis","Cutibacterium_granulosum","skin_imperfection_Para40","Face_Pigmentation_Frequency","Nicotine","Cotinine","Acenaphtylene","Acenaphtene","Fluorene",
                                                 "Phenanthrene","Anthracene","Fluoranthene","Pyrene","Benzo_a_anthracene","Chrysene","Benzo_b_fluoranthene","Benzo_k_fluoranthene","Benzo_a_pyrene","Indeno_1_2_3_cd_pyrene","Dibenzo_a_h_anthracene","Benzo_ghi_perylene")]
Meta_2 <- read.tidy("PCoA_coordinates_kraken2_w_ctypes.txt")[c("Sample","Type")]
Meta <- merge(Meta,Meta_2,by="Sample",all.X=TRUE)
Meta <- Meta[match(labels(BrayCurtis), Meta$Sample), ]

dbRDA <- capscale(BrayMatrix ~ Nicotine+Cotinine+Acenaphtylene+Acenaphtene+Fluorene+
                  Phenanthrene+Anthracene+Fluoranthene+Pyrene+Benzo_a_anthracene+Chrysene+Benzo_b_fluoranthene+Benzo_k_fluoranthene+Benzo_a_pyrene+Indeno_1_2_3_cd_pyrene+Dibenzo_a_h_anthracene+Benzo_ghi_perylene, Meta)
Biplot <- plot(dbRDA)

SampleCoords <- as.data.frame(Biplot$sites)
SampleCoords$Sample <- factor(rownames(SampleCoords))
rownames(SampleCoords) <- NULL
SampleCoords <- merge(SampleCoords, Meta)

EnvCoords <- as.data.frame(Biplot$biplot) * 10
EnvCoords$Variable <- rownames(EnvCoords)
rownames(EnvCoords) <- NULL


Plot <- ggplot(SampleCoords, aes(x = CAP1, y = CAP2))
Plot <- Plot + geom_point(size=2,aes(colour=Type))
Plot <- Plot + theme_bw()
Plot <- Plot + geom_segment(data = EnvCoords, xend = 0, yend = 0, alpha = 0.5, linetype="dotted")
#Plot <- Plot + geom_point(data = EnvCoords, shape=4, size=3)
Plot <- Plot + theme(legend.title = element_text(size = 20))
Plot <- Plot + theme(axis.text = element_text(size = 18))
Plot <- Plot + theme(axis.title = element_text(size = 18, face = "bold"))
Plot <- Plot + theme(legend.text = element_text(size = 14, face = "bold"))
Plot <- Plot + theme(legend.title = element_text(size = 20, face = "bold"))
#Plot <- Plot + scale_colour_brewer(palette="Dark2")
Plot <- Plot + theme(legend.key = element_blank())
Plot <- Plot + theme_classic()
#Optional (add text to lines in dbRDA)
Plot <- Plot + geom_text(size=2,data = EnvCoords, aes(label = Variable))
ggsave("dbRDA_PAH_cutotype_PAM.pdf",width=7,height=5,units="in")

#Compare pairwise distance between and within C-types
Table <- read.tidy("BC_tax_matrix.txt")
Table <- melt(Table)
colnames(Table)[1] <- "Sample1"
colnames(Table)[2] <- "Sample2"
colnames(Table)[3] <- "Bray"

Meta <- read.tidy("PCoA_coordinates_kraken2_w_ctypes.txt")[c("Sample","Type")]
Meta_2 <- read.tidy("../meta/meta_pahoh_new.txt") [c("Sample","City","skin_imperfection_Para40")]
Meta <- merge(Meta,Meta_2,by="Sample",all.X=TRUE)
Merge <- merge(Table, Meta, by.x = "Sample1", by.y = "Sample", all = TRUE) 
names(Merge)[4:6] <- c("Type1","City1","Acne1")
Merge <- merge(Merge, Meta, by.x = "Sample2", by.y = "Sample", all = TRUE)
names(Merge)[7:9] <- c("Type2","City2","Acne2")
#Remove self/self comparisons
Merge <- Merge[which( ! Merge$Sample1 == Merge$Sample2), ]

Merge$Cutotype <- ifelse(Merge$Type1 == Merge$Type2, "Same Type", "Different Type")
Merge$CityType <- ifelse(Merge$City1 == Merge$City2, "Same City","Different City")
Merge$TypeType <- paste0(Merge$Type1, " vs. ", Merge$Type2)
Merge$AcneType <- paste0(Merge$Acne1, " vs. ", Merge$Acne2)

write.tidy(Merge,"pairwise_BC_type.txt")

#Test statistical significance between distances
kruskal.test(Merge$TypeType, Merge$Bray,data=Merge)
Kruskal-Wallis rank sum test

data:  Merge$TypeType and Merge$Bray
Kruskal-Wallis chi-squared = 13813, df = 7625, p-value < 2.2e-16

#Plot BC dissim across cutotypes
Table <- read.tidy("pairwise_BC_type.txt")
Table$TypeType <- factor(Table$TypeType, levels=c("Cutotype 1 vs. Cutotype 1","Cutotype 2 vs. Cutotype 2","Cutotype 2 vs. Cutotype 1"))
Plot <- ggplot(Table,aes(x=TypeType,y=Bray,color=TypeType))
Plot <- Plot + geom_boxplot(outlier.shape=NA) + theme_classic()
ggsave("pairwise_cutotype.pdf",width=4,height=4,units="in")

#Calculate statistical analysis for inter-sample difference significance
library(pgirmess)
kruskal.test(Bray~TypeType,data=Table)
#Kruskal-Wallis rank sum test

#data:  Bray by TypeType
#Kruskal-Wallis chi-squared = 7000.4, df = 2, p-value < 2.2e-16

kruskalmc(Bray~TypeType,data=Table)
Multiple comparison test after Kruskal-Wallis 
#p.value: 0.05 
#Comparisons
#obs.dif critical.dif difference
#Cutotype 1 vs. Cutotype 1-Cutotype 2 vs. Cutotype 2 4714.0675     192.2901       TRUE
#Cutotype 1 vs. Cutotype 1-Cutotype 2 vs. Cutotype 1 5493.0316     170.6790       TRUE
#Cutotype 2 vs. Cutotype 2-Cutotype 2 vs. Cutotype 1  778.9641     204.8148       TRUE