
#
#Sanjeev Sariya
#06 14 2021

.libPaths(c( "/home/ss5505/libraries_R/R_LIB4.0",.libPaths()))

library(dplyr)
library(ggplot2)
library(ggpubr)

df_zf_human_orthologues<-read.table("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/model_organisms/zebra_fish/DEG_files_ZF/human2zebrafish.txt", header=TRUE)

load(file = "/mnt/mfs/hgrcgrid/shared/GT_ADMIX/model_organisms/zebra_fish/DEG_files_ZF/dre_AB42vsCTL.rda")

dre_AB42vsCTL$ENS_id<-rownames(dre_AB42vsCTL)

df_mayo_deg<- read.table("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/model_organisms/zebra_fish/DEG_files_human/STable4_MAYO_SumStats.csv", header=TRUE,sep=",")
df_ch_deg<- read.table("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/model_organisms/zebra_fish/DEG_files_human/STable1_CH_SumStats.csv", header=TRUE,sep=",")
df_rosmap_deg<- read.table("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/model_organisms/zebra_fish/DEG_files_human/STable3_ROSMAP_SumStats.csv", header=TRUE, sep=",")

colnames(df_mayo_deg)<- paste0("mayo_",colnames(df_mayo_deg))
colnames(df_rosmap_deg)<- paste0("rosmap_",colnames(df_rosmap_deg))
colnames(df_ch_deg)<- paste0("CH_",colnames(df_ch_deg))

merge_tstats<-function(temp_df1, temp_df2, col1){
orthol_merged<- merge( temp_df1, df_zf_human_orthologues, by.x = col1, by.y="Gene_humanstable_ID"  )
return ( merge(orthol_merged, dre_AB42vsCTL , by.x ="Gene_ZFstable_ID", by.y="ENS_id"))
}


##merge using the functions
ch_mergeded<-merge_tstats(df_ch_deg, dre_AB42vsCTL,"CH_gene" )
mayo_mergeded<-merge_tstats(df_mayo_deg, dre_AB42vsCTL,"mayo_gene" )
rosmap_mergeded<-merge_tstats(df_rosmap_deg, dre_AB42vsCTL,"rosmap_gene" )

mayo_mergeded_filtered<-mayo_mergeded[which( !is.na(mayo_mergeded$stat) &  !is.na(mayo_mergeded$mayo_t )),]
rosmap_mergeded_filtered<-rosmap_mergeded[which( !is.na(rosmap_mergeded$stat) &  !is.na(rosmap_mergeded$rosmap_t )),]
ch_mergeded_filtered<-ch_mergeded[which( !is.na(ch_mergeded$stat) &  !is.na(ch_mergeded$CH_t)),]

ch_mergeded_filtered_pvalue<-ch_mergeded_filtered[which(ch_mergeded_filtered$CH_P.Value <= 0.05),]
rosmap_mergeded_filtered_pvalue <-rosmap_mergeded_filtered[which(rosmap_mergeded_filtered$rosmap_P.Value  <= 0.05), ]
mayo_mergeded_filtered_pvalue<-mayo_mergeded_filtered[which(mayo_mergeded_filtered$mayo_P.Value <=0.05),]

cor.test( ch_mergeded_filtered_pvalue$CH_t, ch_mergeded_filtered_pvalue$stat, method="spearman")
cor.test( mayo_mergeded_filtered_pvalue$mayo_t, mayo_mergeded_filtered_pvalue$stat, method="spearman")
cor.test( rosmap_mergeded_filtered_pvalue$rosmap_t, rosmap_mergeded_filtered_pvalue$stat, method="spearman")

ch_mergeded_filtered_pvalue_zf<-ch_mergeded_filtered[which(ch_mergeded_filtered$CH_P.Value <= 0.05 & ch_mergeded_filtered$pvalue <=0.05),]
rosmap_mergeded_filtered_pvalue_zf <-rosmap_mergeded_filtered[which(rosmap_mergeded_filtered$rosmap_P.Value  <= 0.05 & rosmap_mergeded_filtered$pvalue<=0.05), ]
mayo_mergeded_filtered_pvalue_zf<-mayo_mergeded_filtered[which(mayo_mergeded_filtered$mayo_P.Value <=0.05 & mayo_mergeded_filtered$pvalue<=0.05),]

cor.test( ch_mergeded_filtered_pvalue_zf$CH_t, ch_mergeded_filtered_pvalue_zf$stat, method="spearman")
cor.test( mayo_mergeded_filtered_pvalue_zf$mayo_t, mayo_mergeded_filtered_pvalue_zf$stat, method="spearman")
cor.test( rosmap_mergeded_filtered_pvalue_zf$rosmap_t, rosmap_mergeded_filtered_pvalue_zf$stat, method="spearman")


bitmap("CH_zf_rnaseq_overlappingsign_genes.tiff")
ggscatter(ch_mergeded_filtered_pvalue_zf, x = "CH_t", y = "stat" ) +  stat_cor(method = "spearman") 
dev.off()

bitmap("mayo_zf_rnaseq_overlappingsign_genes.tiff")
ggscatter(mayo_mergeded_filtered_pvalue_zf, x = "mayo_t", y = "stat" ) +  stat_cor(method = "spearman") 
dev.off()

bitmap("rosmap_zf_rnaseq_overlappingsign_genes.tiff")
ggscatter(rosmap_mergeded_filtered_pvalue_zf, x = "rosmap_t", y = "stat" ) +  stat_cor(method = "spearman") 
dev.off()
