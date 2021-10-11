#Open significant maaslin2 results

Table <- read.tidy("significant_results_for_maaslin_humann_pathway_plot.txt")

Plot <- ggplot(Table,aes(x = reorder(feature, -coef),y=coef,fill=value))
Plot <- Plot + geom_bar(stat="identity") + scale_y_continuous(expand=c(0,0)) + coord_flip() + theme_classic()

ggsave("maaslin2_output_cutotype.pdf",width=5,height=6,units="in")
