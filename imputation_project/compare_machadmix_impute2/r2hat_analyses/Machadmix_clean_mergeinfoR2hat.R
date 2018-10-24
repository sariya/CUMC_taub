#!/bin/Rscript

#10/24/2018
#Sanjeev Sariya
#clean R2hat and merg MachAdmix Info files
##This is for 1000G reference panel mainly

#need info merged file and R2 hat 
#R2hat tool version r2_hat_v109
#
#Imputation is done by making pat and Mat ids as zero by communication from Author Yun Li - date 10/23/2018
#

library(dplyr)
#cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/accuracy_imputation/r2_hat_calculation/machadmix/calculate_r2hat_perchunk/10242018_calculation

info.file<-"mergedInfo_CHR21_MachAdmix"
file.r2.hat<-"merged_r2_hat_MachAdmixCHR21"

df.r2hat<-data.table::fread(file.r2.hat, showProgress = TRUE,header=TRUE)
df.info<-data.table::fread(info.file, showProgress = TRUE,header=TRUE)
print(dim(df.r2hat))
print(dim(df.info))

#--find snps with - value in df.r2hat
print(length(which(df.r2hat$R2_hat=="-")))

print(df.r2hat[head(which(df.r2hat$R2_hat=="-")),])

#- get rid of snps that have - as r2hat 
df.r2hat_cleaned<-df.r2hat[-c(which(df.r2hat$R2_hat=="-")),]
print(dim(df.r2hat_cleaned))

Al1_length<-unlist(lapply(df.info$Al1, function(x) nchar(x)))
Al2_length<-unlist(lapply(df.info$Al2, function(x) nchar(x)))

###throw whatver was an structual variants
df.info<-df.info[which(Al1_length==1 & Al2_length==1),]
print(dim(df.info)) #  598943 SNPs

########
joined_info_r2hat<-left_join(df.info,df.r2hat_cleaned,by=c("SNP"="Marker"))
print(dim(joined_info_r2hat))
print(dim(joined_info_r2hat[complete.cases(joined_info_r2hat),]))

joined_info_r2hat<-joined_info_r2hat[complete.cases(joined_info_r2hat),]
print(dim(joined_info_r2hat)) #584323

joined_info_r2hat$R2_hat<-as.numeric(joined_info_r2hat$R2_hat)
print(min(joined_info_r2hat$R2_hat))
print(max(joined_info_r2hat$R2_hat))

write.table(joined_info_r2hat,"joined_r2hat_machadmixinfo_CHR21_cleaned",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

###
####
"

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
 [1] tidyselect_0.2.5  compiler_3.4.2    magrittr_1.5      assertthat_0.2.0
 [5] R6_2.3.0          pillar_1.3.0      bindrcpp_0.2.2    glue_1.3.0
 [9] tibble_1.4.2      crayon_1.3.4      Rcpp_0.12.19      data.table_1.11.8
[13] pkgconfig_2.0.2   rlang_0.3.0       purrr_0.2.5       bindr_0.1.1

"