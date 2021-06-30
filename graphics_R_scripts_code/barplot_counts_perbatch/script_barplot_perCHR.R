
##
##Sanjeev Sariya
##06 30 2021

.libPaths(c( "C:/Users/ss5505/Documents/R_learning/manual_libraries",.libPaths()))
library(ggplot2)
library(resahpe2)
library(grid)

counts_df<-read.table( "counts_perCHR_perbatch.txt", header=TRUE)
counts_df$CHR<- as.factor(counts_df$CHR)
 
molten_chr_counts<-reshape2::melt(counts_df)

##we need chromosome in text strip
molten_chr_counts$CHROM <-paste("Chromosome",(molten_chr_counts$CHR))

##toughest part in this plot is to get the order
##https://stackoverflow.com/a/11282791/2740831
molten_chr_counts$CHROM_reorder <- factor(molten_chr_counts$CHROM, levels = unique(molten_chr_counts$CHROM[order(molten_chr_counts$CHR)]))

molten_chr_counts$value_per_k <-molten_chr_counts$value/1000
 #  ggplot(molten_chr_counts, aes(fill=variable, y=value, x=CHR)) +  ggplot2::geom_bar(position="dodge", stat="identity")

plot_gg<- ggplot(molten_chr_counts , aes(x = variable, y = value_per_k, fill=as.factor(variable))) +
labs(fill = "Batch") +   ggplot2::geom_bar(position="dodge", stat="identity") +
facet_wrap(~CHROM_reorder,nrow=5 ,scales = "free") +
# scale_y_continuous(limits=c(0,100),breaks=c(10,35,55,100))  +
labs(title="" ,x ="", y = "HQ SNP Count (in thousands)") +     
theme_bw() + theme(plot.title = element_text( size=15, face="bold", hjust = 0.5),
axis.text.x=element_blank(), axis.ticks.x=element_blank(),
legend.background = element_rect(size=0.5,linetype="solid", colour ="black"),  ## , colour ="black"
legend.text=element_text(size=14),legend.title=element_text(size=16),
strip.text.x = element_text(size = 10, face="bold"),
axis.title.x = element_text(size=20,face="bold"),
axis.title.y = element_text(size=20,face="bold") ,
strip.background = element_rect( size = 1.1),panel.border=element_rect(size=1.2)
)

####################exit
####################exit
####################exit
####################exit

#plot_tab <- ggplotGrob(plot_gg)

## https://stackoverflow.com/a/52028670
# gtable_filter_remove <- function (x, name, trim = TRUE){
    matches <- !(x$layout$name %in% name)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    if (trim)
        x <- gtable_trim(x)
    x
}
## Functino ends to remove labels 

#p_filtered <- gtable_filter_remove(plot_tab , name = paste0("axis-l-", c(2, 3,4,5), "-2"),   trim = FALSE)
#p_filtered <- gtable_filter_remove(p_filtered , name = paste0("axis-l-", c(2, 3,4,5), "-3"), trim = FALSE)
#p_filtered <- gtable_filter_remove(p_filtered, name = paste0("axis-l-", c(2, 3,4,5), "-4"),  trim = FALSE)
#p_filtered <- gtable_filter_remove(p_filtered, name = paste0("axis-l-", c(2, 3,4,5), "-5"),  trim = FALSE)
#p_filtered <- gtable_filter_remove(p_filtered, name = paste0("axis-l-", c(1), "-5"), trim = FALSE)

#p_filtered <- gtable_filter_remove(p_filtered, name = paste0("axis-l-", c(1), "-2"),  trim = FALSE)
#p_filtered <- gtable_filter_remove(p_filtered, name = paste0("axis-l-", c(1), "-3"),  trim = FALSE)
#p_filtered <- gtable_filter_remove(p_filtered, name = paste0("axis-l-", c(1), "-4"),  trim = FALSE)

#CairoTIFF(filename = "averageLA_perCHR_acrossGenome.tiff", height=800, width=900)

#grid.newpage()
#grid.draw(p_filtered)
#dev.off()


