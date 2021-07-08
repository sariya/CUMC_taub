
.libPaths(c( "/home/ss5505/libraries_R/R_LIB4.0",.libPaths()))
library(data.table)
library(Cairo)
library(ggplot2)

prepare_dt<-function(file_path){

chr_list<-list()
for(chr in seq(1:22) ) {
    file<-( paste0(file_path,"/", "CHR",chr,"_genotyped_snps") )
    chr_list[[chr]] <- as.data.table(lapply(file, data.table::fread))
chr_list[[chr]][,c("CHR") := chr] 
}
dt_merged <- rbindlist( chr_list)
dt_merged[,c("ALT_AF_diff_MAF") := ALT_Frq - MAF] 

dt_merged$CHROM <-paste("Chromosome",(dt_merged$CHR))

dt_merged$CHROM_reorder <- factor(dt_merged$CHROM, levels = unique(dt_merged$CHROM[order(dt_merged$CHR)]))

return(dt_merged)

}

illumina_flipped_dt<-prepare_dt("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/TOPMED_imputation_analyses/imputation_input_files_perbatch/BATCH2/topmed_output/perCHR_genotype_info_files")

colnames(illumina_flipped_dt) <-c( "SNP"  ,"REF.0.","ALT.1." , "ALT_Frq_illumina_flipped",
"MAF_illumina_flipped" ,"AvgCall_illumina_flipped","Rsq_illumina_flipped",
"Genotyped","LooRsq_illumina_flipped", "EmpR_illumina_flipped", 
"EmpRsq_illumina_flipped","Dose0_illumina_flipped","Dose1_illumina_flipped","CHR", "ALT_AF_diff_MAF_illumina_flipped" ,"CHROM","CHROM_reorder")

illumina_nonflipped_dt<-prepare_dt("/mnt/mfs/hgrcgrid/shared/HISPANICS_IMPUTATION_2017/TOPMED/BATCH2/TOPMED_final_output/perCHR_genotype_info_files/")

colnames(illumina_nonflipped_dt) <-c( "SNP"  ,"REF.0.","ALT.1." , "ALT_Frq_illumina_nonflipped",
"MAF_illumina_nonflipped" ,"AvgCall_illumina_nonflipped","Rsq_illumina_nonflipped",
"Genotyped","LooRsq_illumina_nonflipped", "EmpR_illumina_nonflipped", 
"EmpRsq_illumina_nonflipped","Dose0_illumina_nonflipped","Dose1_illumina_nonflipped",
"CHR", "ALT_AF_diff_MAF_illumina_nonflipped" ,"CHROM","CHROM_reorder")

illumina_nonflipped_dt[,c("CHR","CHROM_reorder","CHROM","Genotyped") :=NULL]

length(intersect( illumina_nonflipped_dt$SNP , illumina_flipped_dt$SNP))
merged_overlappingSNPs<-illumina_nonflipped_dt[illumina_flipped_dt, on  = "SNP", nomatch=0]
dim(merged_overlappingSNPs)

text_settings<-theme_bw() + theme( plot.title = element_text( size=15, face="bold", hjust = 0.5),
strip.text.x = element_text(size = 10, face="bold"),
axis.title.x = element_text(size=20,face="bold"),
axis.title.y = element_text(size=20,face="bold") 
)

##plot rsq 

plot_rsq<- ggplot2::ggplot(merged_overlappingSNPs,aes(x = Rsq_illumina_flipped, y =Rsq_illumina_nonflipped)) +  
geom_point(size=1) + facet_wrap(~CHROM_reorder,nrow=5 ,scales = "free") + labs(title="" ,x ="Rsq illumina_flipped", y = "Rsq illumina_nonflipped") + text_settings

CairoTIFF(filename = "rsq_perCHR_comp.tiff", height=800, width=900)
plot_rsq
dev.off()

##plot diff ALT AF and MAF

plot_diff_alt_MAF<- ggplot2::ggplot(merged_overlappingSNPs,aes(x = Rsq_illumina_flipped, y =Rsq_illumina_nonflipped)) +  
geom_point(size=1) + facet_wrap(~CHROM_reorder,nrow=5 ,scales = "free") + labs(title="" ,x ="Diff ALT-MAF illumina_flipped", y = "Diff ALT-MAF illumina_nonflipped") + text_settings

CairoTIFF(filename = "diff_alt_MAF_perCHR_comp.tiff", height=800, width=900)
plot_diff_alt_MAF
dev.off()

##plot emp RSQ

plot_empRSQ<- ggplot2::ggplot(merged_overlappingSNPs,aes(x = EmpRsq_illumina_flipped, y =EmpRsq_illumina_nonflipped)) +  
geom_point(size=1) + facet_wrap(~CHROM_reorder,nrow=5 ,scales = "free") + labs(title="" ,x ="EmpRsq illumina_flipped", y = "EmpRsq illumina_nonflipped") + text_settings

CairoTIFF(filename = "empRSQ_perCHR_comp.tiff", height=800, width=900)
plot_empRSQ
dev.off()

