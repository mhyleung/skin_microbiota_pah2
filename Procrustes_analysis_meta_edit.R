library(wilkoxmisc)
library(ggplot2)
library(devtools)
library(vegan)
library(phyloseq)
library(Rmisc)
library(readr)
library(dplyr)
library(RVAideMemoire)
library(ade4)


#Kraken2 vs. Humann3
tax_cast <- read.table("bracken2_all_samples_wide_0_00001_adjusted_for_human_clean_no_con_122.txt",header=T,row.names=1)
tax_cast <- t(tax_cast)
distance_b_kraken <- vegdist(tax_cast, method="bray")

func <- read.table("joined_genefamilies_table_unstratified_cpm_ko_122.txt" ,header=T,row.names=1)
tax_cast_function <- t(func)
distance_b_function <- vegdist(tax_cast_function, method="bray")

mantel <- mantel.rtest(distance_b_kraken, distance_b_function, nrepet = 999)
#Monte-Carlo test
#Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)

#Observation: 0.6855185 

#Based on 999 replicates
#Simulated p-value: 0.001 
#Alternative hypothesis: greater 

#Std.Obs Expectation    Variance 
#9.962840566 0.001977615 0.004707200
procrustes <- protest(distance_b_kraken, distance_b_function, nperm = 999)
#Call:
#  protest(X = tax_BC, Y = function_BC, nperm = 999) 

#Procrustes Sum of Squares (m12 squared):        0.4441 
#Correlation in a symmetric Procrustes rotation: 0.7456 
#Significance:  0.001 

#Permutation: free
#Number of permutations: 999


data1 = data.frame(rda1 = pro.test$Yrot[,1], rda2 = pro.test$Yrot[,2], 
                   xrda1 = pro.test$X[,1], xrda2 = pro.test$X[,2])
write_tsv(data1, "data1_tax_BC.tsv")
data2 = data.frame(rda1 =  pro.test$X[,1], rda2 = pro.test$X[,2], 
                   xrda1 = pro.test$Yrot[,1], xrda2 = pro.test$Yrot[,2])
write_tsv(data2,"data2_function_BC.tsv")

data1 <- merge(data1, meta, by = "Sample")
data2 <- merge(data2, meta, by = "Sample")


data1 <- read_tsv("data1_tax_BC.tsv")
data2 <- read_tsv("data2_function_BC.tsv")
data = rbind(data1, data2)
##data$type = c(c(rep("ASV", 131), rep("OTU", 131)),
#c(rep("ASV", 131), rep("OTU", 131)))




pro.test_1 = procrustes(distance_b_kraken, distance_b_function, scale = TRUE)
dataA1 = data.frame(rda1 = pro.test_1$Yrot[,1], rda2 = pro.test_1$Yrot[,2], 
                   xrda1 = pro.test_1$X[,1], xrda2 = pro.test_1$X[,2])
write_tsv(dataA1, "data1_tax_bc.tsv")
dataA2 = data.frame(rda1 =  pro.test_1$X[,1], rda2 = pro.test_1$X[,2], 
                   xrda1 = pro.test_1$Yrot[,1], xrda2 = pro.test_1$Yrot[,2])
write_tsv(dataA2,"data2_function_bc.tsv")

dataA1 <- read.tidy("data1_tax_BC.tsv")
dataA2 <- read.tidy("data2_function_BC.tsv")
meta <- read_tsv("../beta_adonis2//PCoA_coordinates_kraken2_w_ctypes.txt")[c("Sample","Type")]
dataA1 <- merge(dataA1, meta, by = "Sample", all.x=TRUE)
dataA2 <- merge(dataA2, meta, by = "Sample", all.x=TRUE)
dataA = rbind(dataA1, dataA2)

cutoff_x = data.frame( x = c(-Inf, Inf), y = 0, cutoff_x = factor(50) )
cutoff_y = data.frame( y = c(-Inf, Inf), x = 0, cutoff_y = factor(50) )
PlotA <- ggplot(dataA) +
  geom_line(aes( x, y, linetype = cutoff_x ), cutoff_x, show.legend = FALSE) +
  geom_line(aes( x, y, linetype = cutoff_y ), cutoff_y, show.legend = FALSE)

PlotA <- PlotA +
  geom_point(aes(x = rda1,  y = rda2,  shape = Method, colour = Type), size = 3, alpha=0.7)

PlotA <- PlotA + scale_color_manual(values = c("blue", "red"))

PlotA <- PlotA + geom_segment(aes(x = rda1,y = rda2, xend = xrda1,yend = xrda2), size=0.2,alpha=0.5)
PlotA <- PlotA + theme_bw() + theme(legend.position = "right") + xlab("Dimension 1") + ylab("Dimension 2") +
  theme(axis.ticks = element_blank(),
        axis.title = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))
ggsave("procrustes_kraken2_function_BC_cutotype.pdf", device = "pdf", width = 8, height = 8, units="in")
