#draw SMEG heatmap
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(wilkoxmisc)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(gplots)
library(heatmap3)

#Cacnes
CLong <- read.tidy("c_acnes_smeg_long.txt")
Table <- CLong[which(CLong$SMEG >=1) ,]
Table <- dcast(Table, Cluster ~ Sample, value.var = "SMEG")

rownames(Table) <- Table$Cluster
Table$Cluster <- NULL

Matrix <- as.matrix(Table)
write.tidy(Matrix,"c_acnes_smeg_matrix.txt",row.names=TRUE)


Matrix <- read.table("c_acnes_smeg_matrix_edited_no_rows_cols.txt")
Matrix <- as.matrix(Matrix)
Cluster <- read.tidy("c_acnes_smeg_matrix_cluster.txt")
Cuto <- read.tidy("c_acnes_smeg_matrix_cutotype.txt")

rownames(Matrix) <- Cluster$Cluster
colnames(Matrix) <- Cuto$Cutotype

#Meta <- read.tidy("../meta/meta_pahoh_newer.txt")
#Meta <- as.data.frame(Meta)

my_group <- as.factor(Cuto$Cutotype)
colSide <- brewer.pal(8,"Set3")[my_group]
colMain <- colorRampPalette(c("black", "red","orange","yellow"))(299) #Set main heatmap colour
heatmap3(Matrix,Colv=NA, Rowv=NA, scale="none", margins = c(1, 6),cexRow=0.75,labCol=NA, ColSideColors=colSide, col=colMain)
heatmap3(Matrix, showRowDendro=FALSE, labCol=NA, ColSideColors=colSide, ColSideLabs="Cutotype", col=colMain)


#Save_export

Merge <- merge(Long,Meta,by="Sample",all.x=TRUE)
Cluster4 <- Merge[which(Merge$Cluster == "Cluster 4") ,]
wilcox.test(SMEG~Acne_Status,data=Cluster4)

#Not statistically significant for all clusters

#Mluteus
MLong <- read.tidy("m_luteus_smeg_long.txt")
Table <- MLong[which(MLong$SMEG >=1) ,]
Table <- dcast(Table, Cluster ~ Sample, value.var = "SMEG")

rownames(Table) <- Table$Cluster
Table$Cluster <- NULL
Matrix <- as.matrix(Table)
write.tidy(Matrix,"m_luteus_smeg_matrix.txt",row.names=TRUE)

Matrix <- read.table("m_luteus_smeg_matrix_no_col_no_row.txt")
Matrix <- as.matrix(Matrix)
Cluster <- read.tidy("m_luteus_smeg_matrix_cluster.txt")
Cuto <- read.tidy("m_luteus_smeg_matrix_cutotype.txt")

rownames(Matrix) <- Cluster$Cluster
colnames(Matrix) <- Cuto$Cutotype

#Meta <- read.tidy("../meta/meta_pahoh_newer.txt")
#MLong <- read.tidy("m_luteus_smeg_long.txt")[c("Sample","City")]
#MLong <- unique(MLong)
#Meta <- merge(MLong,Meta,by="Sample",all.X=TRUE)
#Meta <- as.data.frame(Meta)

my_group_new <- as.factor(Meta$City.x)
colSide <- brewer.pal(8,"Set3")[my_group_new]
colMain <- colorRampPalette(c("black", "red","orange","yellow"))(299) #Set main heatmap colour
heatmap3(Matrix,Colv=NA, Rowv=NA, scale="none", margins = c(1, 6),cexRow=0.75,labCol=NA, ColSideColors=colSide, col=colMain)
#Save_export