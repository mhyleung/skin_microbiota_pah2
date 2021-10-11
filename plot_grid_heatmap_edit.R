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

Table <- read_tsv("GRID_input_new_p10.txt")
Table <- melt(Table, id.vars = "Genome", value.name="GRID")


#Add cutotype information
Meta <- read.tidy("../../beta_adonis2/PCoA_coordinates_kraken2_w_ctypes.txt")[c("Sample","Type")]
Table <- merge(Table,Meta,by.x="variable",by.y="Sample",all.x=TRUE)
write.tidy(Table,"GRID_input_new_p10_long.txt")
Table <- read.tidy("GRID_input_new_p10_long.txt")
Table <- Table[which(Table$GRID >=1) ,]
Table <- dcast(Table, Genome ~ variable, value.var = "GRID")

rownames(Table) <- Table$Genome
Table$Genome <- NULL
Matrix <- data.matrix(Table)

write.tidy(Matrix,"test_matrix.txt",row.names=TRUE)

Matrix <- read.table("test_matrix_no_row_col.txt")
Matrix <- as.matrix(Matrix)

Species <- read.tidy("matrix_order_species.txt")
Cuto <- read.tidy("matrix_order_cutotype.txt")

rownames(Matrix) <- Species$Species
colnames(Matrix) <- Cuto$Cutotype

my_group <- as.factor(Cuto$Cutotype)
colSide <- brewer.pal(3,"Set3")[my_group] #add color to factor
colMain <- colorRampPalette(c("black", "red","orange","yellow"))(299) #Set main heatmap colour
heatmap3(Matrix, Colv=NA, Rowv=NA, scale="none", margins = c(1, 11),cexRow=0.75,labCol=NA, ColSideColors=colSide, ColSideLabs="Cutotype", col=colMain)
quartz.save(type="pdf",dpi=1200,file="grid_heatmap_by_cutotype.pdf")
