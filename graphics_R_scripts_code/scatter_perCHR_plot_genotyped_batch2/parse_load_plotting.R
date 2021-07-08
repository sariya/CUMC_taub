
##
##07 08 2021
##Sanjeev Sariya

.libPaths(c( "/home/ss5505/libraries_R/R_LIB4.0",.libPaths()))
library(data.table)
library(Cairo)
library(ggplot2)

chr_list<-list()

for(chr in seq(1:22) ) {
    file<-( paste0("CHR",chr,"_genotyped_snps") )
    chr_list[[chr]] <- as.data.table(lapply(file, data.table::fread))
chr_list[[chr]][,c("CHR") := chr] 
}
dt_merged <- rbindlist( chr_list)
dt_merged[,c("ALT_AF_diff_MAF") := ALT_Frq - MAF] 

dt_merged$CHROM <-paste("Chromosome",(dt_merged$CHR))

##set text - font settings
text_settings<-theme_bw() + theme( plot.title = element_text( size=15, face="bold", hjust = 0.5),
strip.text.x = element_text(size = 10, face="bold"),
axis.title.x = element_text(size=20,face="bold"),
axis.title.y = element_text(size=20,face="bold") 
)

dt_merged$CHROM_reorder <- factor(dt_merged$CHROM, levels = unique(dt_merged$CHROM[order(dt_merged$CHR)]))
plot_gg<- ggplot2::ggplot(dt_merged,aes(x = Rsq)) + geom_histogram(position="identity", alpha=0.5) +
facet_wrap(~CHROM_reorder,nrow=5 ,scales = "free") + labs(title="" ,x ="Rsq imputed", y = "") + text_settings

#theme_bw() + theme( plot.title = element_text( size=15, face="bold", hjust = 0.5), #strip.text.x = element_text(size = 10, face="bold"),
#axis.title.x = element_text(size=20,face="bold"), #axis.title.y = element_text(size=20,face="bold")  #)

CairoTIFF(filename = "rsq_perCHR.tiff", height=800, width=900)
plot_gg
dev.off()

plot_gg_rsq_alt<- ggplot2::ggplot(dt_merged,aes(x = Rsq, y =ALT_Frq)) +  geom_point(size=1) +
facet_wrap(~CHROM_reorder,nrow=5 ,scales = "free") + labs(title="" ,x ="Rsq imputed", y = "ALT_Frq") + text_settings

CairoTIFF(filename = "rsqALT_freq_perCHR.tiff", height=800, width=900)
plot_gg_rsq_alt
dev.off()

#axis.text.x=element_blank(),
#axis.ticks.x=element_blank(),

plot_gg_rsq_alt_loo<- ggplot2::ggplot(dt_merged,aes(x = Rsq, y =LooRsq)) +  geom_point(size=1) +
facet_wrap(~CHROM_reorder,nrow=5 ,scales = "free") + labs(title="" ,x ="Rsq imputed", y = "Leave one out Emp Rsq") + text_settings

CairoTIFF(filename = "rsq_looRSQ_perCHR.tiff", height=800, width=900)
plot_gg_rsq_alt_loo
dev.off()

gg_altfreq_alt_loo<- ggplot2::ggplot(dt_merged,aes(x = ALT_Frq, y =LooRsq)) +  
geom_point(size=1) + facet_wrap(~CHROM_reorder,nrow=5 ,scales = "free") + labs(title="" ,x ="ALT Frq", y = "Leave one out Emp Rsq") + text_settings

CairoTIFF(filename = "altFreq_looRSQ_perCHR.tiff", height=800, width=900)
gg_altfreq_alt_loo
dev.off()

gg_altfreq_empRSQ<- ggplot2::ggplot(dt_merged,aes(x = ALT_Frq, y =EmpRsq)) +  
geom_point(size=1) + facet_wrap(~CHROM_reorder,nrow=5 ,scales = "free") + labs(title="" ,x ="ALT Frq", y = "EmpRsq") + text_settings

CairoTIFF(filename = "altFreq_empRSQ_perCHR.tiff", height=800, width=900)
gg_altfreq_empRSQ

dev.off()


plot_AF_diff_emprsq<- ggplot2::ggplot(dt_merged,aes(x = ALT_AF_diff_MAF , y= EmpRsq)) + 
geom_point(size=1) +
facet_wrap(~CHROM_reorder,nrow=5 ,scales = "free") + labs(title="" ,x ="Diff of alt AF with MAF", y = "EmpRsq") + text_settings

CairoTIFF(filename = "alt_diff_maf_empRSQ_perCHR.tiff", height=800, width=900)
plot_AF_diff_emprsq
dev.off()

