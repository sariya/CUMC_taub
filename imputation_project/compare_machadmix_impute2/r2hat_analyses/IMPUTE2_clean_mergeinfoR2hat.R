#!/bin/Rscript

#
#Sanjeev Sariya 10/24/2018
#R2 hat
#This is for 1000G reference panel mainly
#need info merged file and R2 hat 
#R2hat tool version r2_hat_v109

library(dplyr)
info.file<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/impute2/1KGP/merged_files/unharmonized/info/merged_info_CHR21.info"
file.r2.hat<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/accuracy_imputation/r2_hat_calculation/impute2/1000G/calculate_r2hat_perchunk/final_probs/imputed_r2hat.out.r2_hat"

df.r2hat<-data.table::fread(file.r2.hat, showProgress = TRUE,header=TRUE)
df.info<-data.table::fread(info.file, showProgress = TRUE,header=TRUE)
print(dim(df.r2hat))
print(dim(df.info))

#--find snps with - value in df.r2hat
print(length(which(df.r2hat$R2_hat=="-")))
print(head(which(df.r2hat$R2_hat=="-")))
print(df.r2hat[head(which(df.r2hat$R2_hat=="-")),])

#- get rid of snps that have - as r2hat 
df.r2hat_cleaned<-df.r2hat[-c(which(df.r2hat$R2_hat=="-")),]

#--fix MAF for Info

indices_dfref<-df.info$exp_freq_a1>0.5
df.info$exp_freq_a1[indices_dfref]<-(1-df.info$exp_freq_a1)[indices_dfref]

#--get rid of indels from info file

a0_length<-unlist(lapply(df.info$a0, function(x) nchar(x)))
a1_length<-unlist(lapply(df.info$a1, function(x) nchar(x)))

###throw whatver was an structual variants
df.info<-df.info[which(a1_length==1 & a0_length==1),]
print(dim(df.info)) # 1055422 SNPs

#
joined_info_r2hat<-left_join(df.info,df.r2hat_cleaned,by=c("rs_id"="Marker"))
print(dim(joined_info_r2hat))
print(dim(joined_info_r2hat[complete.cases(joined_info_r2hat),]))

joined_info_r2hat<-joined_info_r2hat[complete.cases(joined_info_r2hat),]
print(dim(joined_info_r2hat))

joined_info_r2hat$R2_hat<-as.numeric(joined_info_r2hat$R2_hat)
print(min(joined_info_r2hat$R2_hat))
print(max(joined_info_r2hat$R2_hat))

write.table(joined_info_r2hat,"joined_r2hat_imputeinfo_CHR21_cleaned",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

#df.r2hat_cleaned<-df.r2hat[-c(which(df.r2hat$R2_hat=="-")),]

"
sessionInfo()
sessionInfo()
R version 3.4.2 (2017-09-28)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 9 (stretch)

Matrix products: default
BLAS: /mnt/mfs/cluster/bin/R-3.4/lib/libRblas.so
LAPACK: /mnt/mfs/cluster/bin/R-3.4/lib/libRlapack.so

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] dplyr_0.7.7

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.19      crayon_1.3.4      assertthat_0.2.0  R6_2.3.0
 [5] magrittr_1.5      pillar_1.3.0      rlang_0.3.0       data.table_1.11.8
 [9] bindrcpp_0.2.2    tools_3.4.2       glue_1.3.0        purrr_0.2.5
[13] compiler_3.4.2    pkgconfig_2.0.2   bindr_0.1.1       tidyselect_0.2.5
[17] tibble_1.4.2

"