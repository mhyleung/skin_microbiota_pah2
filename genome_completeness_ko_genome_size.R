#Statistics for correlation between number of KOs and genome size and completeness
library(wilkoxmisc)
Table <- read.tidy("number_ko_genome_completeness.txt")

#All samples
cor.test(Table$KO_No, Table$Completeness, data=Table, method="spearman")
#Spearman's rank correlation rho

#data:  Table$KO_No and Table$Completeness
#S = 2149820, p-value = 0.01499
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1543213 

#Only select hqMAGs
HQ <- Table[which(Table$HQ == "TRUE") ,]

#Spearman ΚΟ vs. genome completeness
cor.test(HQ$KO_No, HQ$Completeness, data=HQ, method="spearman")
#Spearman's rank correlation rho

#data:  HQ$KO_No and HQ$Completeness
#S = 302008, p-value = 0.9823
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#        rho 
#0.002025377 

#Spearman ΚΟ vs. genome size
cor.test(HQ$KO_No, HQ$Size, data=HQ, method="spearman")
#data:  HQ$KO_No and HQ$Size
#S = 54940, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.8184523 