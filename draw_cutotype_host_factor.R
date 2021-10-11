library(wilkoxmisc)
#Draw stacked  bar plot

Table <- read.tidy("cutotype_host_factors.txt")

Plot <- ggplot(Table,aes(x=Variable,y=Percentage)) + geom_col(aes(fill = Cutotype), width = 0.8)
Plot <- Plot + facet_wrap(~Type,scales="free_x") + theme_classic() + scale_y_continuous(expand=c(0,0)) + theme(axis.title.x=element_blank())
Plot <- Plot + ylab(paste0("Prevalence (%)"))
ggsave("cutotype_host.pdf",width=6,height=4,units="in")

#Calculate stats
Age <- Table[which(Table$Type == "Age group") ,]
kruskal.test(Cutotype~Variable,data=Age)
#data:  Cutotype by Variable
#Kruskal-Wallis chi-squared = 0, df = 3, p-value = 1

City <- Table[which(Table$Type == "City") ,]
wilcox.test(Cutotype~Variable,data=City)
#data:  Cutotype by Variable
#Kruskal-Wallis chi-squared = 0, df = 1, p-value = 1

Acne <- Table[which(Table$Type == "Acne") ,]
wilcox.test(Cutotype~Variable,data=City)