#!/bin/Rscript

#
#Date 10/17/2018
#Sanjeev Sariya
#
#Three datasets: A) common HRC-WES, B) 1000G-WES and 
# C) HRC-1000G-WES (for common SNPs in three compare C)a) HRC-WES and C)b) 1000G-WES)
#Once the common are found out: Extract plink data for these SNPs from WES
#make WES as gold standard. 
#Get SNPs with 0% missingness from WES data and extract SNPs for them . make plink --recode A file
#

#For these 100% present common SNPs get geno prob from imputed data for respective reference panel
#Make them into plink with plink --gen geno_prob_100perc --hard-call-threshold 0.1  --make-bed --out plink_dosagebimbam_kgp  --sample kgp_sample.sample
#Make order of WES and the above outputted data same
#Make .raw file using plink --recode A command
#

#Take SNP list present 100% 
#Take .raw from reference panel and WES 
#Get info file to extract/filter SNPs based on AF and Info quality from IMPUTE2

#Need Raw files
#Need info files #need SNPs with 100% presence in Seq data overlapped with resp panel 

library(plyr)
library(dplyr)
options(digits=15)

filter_raw<-function(temp.df,snps_to_keep){
#filter raw file based on input snp list
col_names<-colnames(temp.df)
index<-match(snps_to_keep,col_names)
return (temp.df[,c(1,2,c(index))])
}

####################
##################### Function ends
#####################

filter_info<-function(tempdf_info,temp_snps_forInfo,temp_info){

#get temp df for info, get list of snps for which info is needed and get info threshold
#Return data frame based on info filter

index_snps_matched<-match(temp_snps_forInfo,tempdf_info$snp_harmonized)
tempdf_info<-tempdf_info[index_snps_matched,]
return (tempdf_info[which(tempdf_info$info>=temp_info),10])

}
############Function ends
###
###Function ends

calculate_kappa<-function(dfraw_geno,dfraw_seq) {
#dfraw_geno<-as.data.frame(dfraw_geno) dfraw_seq<-as.data.frame(dfraw_seq)

table_three_bythree<-sapply(seq(3,ncol(dfraw_geno)),function(i) {

#
#Return simplify as FALSE
#
#Create temp matrix to store 3X3 table 
temp.matrix<- matrix(0,3,3)

# Find indices which has 0 1 2
#
index_0<-which(dfraw_seq[,i] ==0)
index_1<-which(dfraw_seq[,i] ==1)
index_2<-which(dfraw_seq[,i] ==2)

#
#Make nine scenarios 
# 0-0,1,2: 1-0,1,2 : 2-0,1,2
##

temp.matrix[1,1]<-sum(dfraw_geno[index_0,i]==0,na.rm =TRUE)
temp.matrix[2,2]<-sum(dfraw_geno[index_1,i]==1,na.rm =TRUE)
temp.matrix[3,3]<-sum(dfraw_geno[index_2,i]==2,na.rm =TRUE)

#
#Ignore NAs in Sum
#

temp.matrix[2,1]<-sum(dfraw_geno[index_0,i]==1,na.rm =TRUE)
temp.matrix[1,2]<-sum(dfraw_geno[index_1,i]==0,na.rm =TRUE)
temp.matrix[1,3]<-sum(dfraw_geno[index_2,i]==0,na.rm =TRUE)

#
#Ignore NAs in Sum
#

temp.matrix[3,1]<-sum(dfraw_geno[index_0,i]==2,na.rm =TRUE)
temp.matrix[3,2]<-sum(dfraw_geno[index_1,i]==2,na.rm =TRUE)
temp.matrix[2,3]<-sum(dfraw_geno[index_2,i]==1,na.rm =TRUE)

return(temp.matrix)
},simplify=F)
####################################################################################

#########################
#####https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
##############################

###Do sum of all matrices
sum_genotype<-Reduce('+', table_three_bythree)

#print(sum_genotype)

if( ( sum(sum_genotype) + sum(is.na(dfraw_geno)) )- ( ncol(dfraw_seq) -2 ) *nrow(dfraw_seq) !=0){
print("We have an issue")
}
###-check difference as good or not


total<-sum(sum_genotype)
agreement_initial_sum<-sum_genotype[1,1]+sum_genotype[2,2]+sum_genotype[3,3]
agreement_initial<-(agreement_initial_sum)/total
print(agreement_initial)

adjust_chance.mat<-matrix(0,3,3)

adjust_chance.mat[1,1]<-(sum(sum_genotype[,1] )*sum(sum_genotype[1,] ))/total
adjust_chance.mat[2,2]<-(sum(sum_genotype[,2] )*sum(sum_genotype[2,] ))/total
adjust_chance.mat[3,3]<-(sum(sum_genotype[,3] )*sum(sum_genotype[3,] ))/total
sum_chance<-sum(adjust_chance.mat)
sum_geno_rows<-sum(apply(sum_genotype,1,sum))

kappa<-(agreement_initial_sum-sum_chance)/(sum_geno_rows-sum_chance)
#print(kappa)

}
############Function ends

###
###Function ends
###

read.raw<-function(temp.file){
##SNPs raw file from plink 
temp.df<-data.table::fread(temp.file, showProgress = TRUE,header=TRUE)
return (as.data.frame(temp.df[,-c(3,4,5,6)]))
}

###
###Function ends
###


#}

###
###Function ends
###

####
infosnps_range<-function(temp.df,info){

print(length(which(temp.df$info>=info)))
return (temp.df[which(temp.df$info>=info),] )

}

###
###Function ends
###

read.info_harmonize<-function(temp.file, chr){

temp_df<-data.table::fread(temp.file, showProgress = TRUE,header=TRUE)
temp_df<-temp_df[,-c(1,9,10)]
temp_df$snp_harmonized<-paste(chr,temp_df$position,temp_df$a0,temp_df$a1,sep=":")

return (temp_df)
}

##############

file.ref.panel.raw<-"100snps_hrc.raw"
file.seq.raw<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/accuracy_imputation/accuracy_WES/correct_keys/last_agreement/100_percent/100seq_seq.raw"
file.info<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/accuracy_imputation/accuracy_WES/correct_keys/merged_hrc_info_CHR14.info"
file.100perc_seqSNPs<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/accuracy_imputation/accuracy_WES/correct_keys/last_agreement/100_percent/100percent"
df.100percsnps<-read.table(file.100perc_seqSNPs,header=FALSE)
df.raw_seq<-read.raw(file.seq.raw)
df.raw_refpanel<-read.raw(file.ref.panel.raw)
df.info_refpanel<-read.info_harmonize(file.info,toString(14))

indices_dfref<-df.info_refpanel$exp_freq_a1>0.5
df.info_refpanel$exp_freq_a1[indices_dfref]<-(1-df.info_refpanel$exp_freq_a1)[indices_dfref]

colnames(df.raw_seq)<- gsub(pattern="_.",replacement="",colnames(df.raw_seq)) 
colnames(df.raw_refpanel)<- gsub(pattern="_.",replacement="",colnames(df.raw_refpanel)) 

print(dim(df.100percsnps))
info_overlapping_snps100perc<-left_join(df.100percsnps,df.info_refpanel,by=c("V1"="snp_harmonized"))
print(dim(info_overlapping_snps100perc))
refp_df_0.4<-infosnps_range(info_overlapping_snps100perc,0.4)
refp_df_0.8<-infosnps_range(info_overlapping_snps100perc,0.8)

refp_df_0_maf_1_5perct<-info_overlapping_snps100perc[which(info_overlapping_snps100perc$exp_freq_a1>=0.01 & info_overlapping_snps100perc$exp_freq_a1<=0.05 ),]
print(dim(refp_df_0_maf_1_5perct))
refp_df_0_maf_.001_1perct<-info_overlapping_snps100perc[which(info_overlapping_snps100perc$exp_freq_a1>=0.001 & info_overlapping_snps100perc$exp_freq_a1< 0.01 ),]
print(dim(refp_df_0_maf_.001_1perct))
refp_df_0_maf_0_.001perct<-info_overlapping_snps100perc[which(info_overlapping_snps100perc$exp_freq_a1>=0 & info_overlapping_snps100perc$exp_freq_a1< 0.001 ),]
print(dim(refp_df_0_maf_0_.001perct))

filtered_raw.ref_0_1_5perc<-filter_raw(df.raw_refpanel,refp_df_0_maf_1_5perct[,1])
dim(filtered_raw.ref_0_1_5perc)
filtered_raw.seq_0_1_5perc<-filter_raw(df.raw_seq,refp_df_0_maf_1_5perct[,1])
dim(filtered_raw.seq_0_1_5perc)

calculate_kappa(filtered_raw.ref_0_1_5perc,filtered_raw.seq_0_1_5perc)

filtered_raw.ref_0_.001_1perct<-filter_raw(df.raw_refpanel,refp_df_0_maf_.001_1perct[,1])
dim(filtered_raw.ref_0_.001_1perct)
filtered_raw.seq_0_.001_1perct<-filter_raw(df.raw_seq,refp_df_0_maf_.001_1perct[,1])
dim(filtered_raw.seq_0_.001_1perct)

calculate_kappa(filtered_raw.ref_0_.001_1perct,filtered_raw.seq_0_.001_1perct)

filtered_raw.ref_0_0_.001perct<-filter_raw(df.raw_refpanel,refp_df_0_maf_0_.001perct[,1])
dim(filtered_raw.ref_0_0_.001perct)
filtered_raw.seq_0_0_.001perct<-filter_raw(df.raw_seq,refp_df_0_maf_0_.001perct[,1])
dim(filtered_raw.seq_0_0_.001perct)

calculate_kappa(filtered_raw.ref_0_0_.001perct,filtered_raw.seq_0_0_.001perct)


###########################
refp_df_0.4_maf_1_5perct<-refp_df_0.4[which(refp_df_0.4$exp_freq_a1>=0.01 & refp_df_0.4$exp_freq_a1<=0.05 ),]
dim(refp_df_0.4_maf_1_5perct)
refp_df_0.4_maf_.001_1perct<-refp_df_0.4[which(refp_df_0.4$exp_freq_a1>=0.001 & refp_df_0.4$exp_freq_a1<0.01 ),]
dim(refp_df_0.4_maf_.001_1perct)
refp_df_0.4_maf_0_.001perct<-refp_df_0.4[which(refp_df_0.4$exp_freq_a1>=0 & refp_df_0.4$exp_freq_a1<0.001 ),]
dim(refp_df_0.4_maf_0_.001perct)

#
#Filter for 80% from 40% DF reference panel
#

refp_df_0.8_maf_1_5perct<-refp_df_0.4_maf_1_5perct[which( refp_df_0.4_maf_1_5perct$info>=0.8),]
dim(refp_df_0.8_maf_1_5perct)
refp_df_0.8_maf_.001_1perct<-refp_df_0.4_maf_.001_1perct[which( refp_df_0.4_maf_.001_1perct$info>=0.8),]
dim(refp_df_0.8_maf_.001_1perct)
refp_df_0.8_maf_0_.001perct<-refp_df_0.4_maf_0_.001perct[which( refp_df_0.4_maf_0_.001perct$info>=0.8),]
dim(refp_df_0.8_maf_0_.001perct)

filtered_raw.ref_0.4_1_5perc<-filter_raw(df.raw_refpanel,refp_df_0.4_maf_1_5perct[,1])
dim(filtered_raw.ref_0.4_1_5perc)
filtered_raw.seq_0.4_1_5perc<-filter_raw(df.raw_seq,refp_df_0.4_maf_1_5perct[,1])
dim(filtered_raw.seq_0.4_1_5perc)

calculate_kappa(filtered_raw.ref_0.4_1_5perc,filtered_raw.seq_0.4_1_5perc)

filtered_raw.ref_0.4_.001_1perct<-filter_raw(df.raw_refpanel,refp_df_0.4_maf_.001_1perct[,1])
dim(filtered_raw.ref_0.4_.001_1perct)
filtered_raw.seq_0.4_.001_1perct<-filter_raw(df.raw_seq,refp_df_0.4_maf_.001_1perct[,1])
dim(filtered_raw.seq_0.4_.001_1perct)

calculate_kappa(filtered_raw.ref_0.4_.001_1perct,filtered_raw.seq_0.4_.001_1perct)

#refp_df_0.4_maf_0_.001perct

filtered_raw.ref_0.4_0_.001perct<-filter_raw(df.raw_refpanel,refp_df_0.4_maf_0_.001perct[,1])
dim(filtered_raw.ref_0.4_0_.001perct)
filtered_raw.seq_0.4_0_.001perct<-filter_raw(df.raw_seq,refp_df_0.4_maf_0_.001perct[,1])
dim(filtered_raw.seq_0.4_0_.001perct)

calculate_kappa(filtered_raw.ref_0.4_0_.001perct,filtered_raw.seq_0.4_0_.001perct)
#########################################################################################

filtered_raw.ref_0.8_1_5perc<-filter_raw(df.raw_refpanel,refp_df_0.8_maf_1_5perct[,1])
dim(filtered_raw.ref_0.8_1_5perc)
filtered_raw.seq_0.8_1_5perc<-filter_raw(df.raw_seq,refp_df_0.8_maf_1_5perct[,1])
dim(filtered_raw.seq_0.8_1_5perc)

calculate_kappa(filtered_raw.ref_0.8_1_5perc,filtered_raw.seq_0.8_1_5perc)

filtered_raw.ref_0.8_.001_1perct<-filter_raw(df.raw_refpanel,refp_df_0.8_maf_.001_1perct[,1])
dim(filtered_raw.ref_0.8_.001_1perct)
filtered_raw.seq_0.8_.001_1perct<-filter_raw(df.raw_seq,refp_df_0.8_maf_.001_1perct[,1])
dim(filtered_raw.seq_0.8_.001_1perct)

calculate_kappa(filtered_raw.ref_0.8_.001_1perct,filtered_raw.seq_0.8_.001_1perct)

#refp_df_0.4_maf_0_.001perct

filtered_raw.ref_0.8_0_.001perct<-filter_raw(df.raw_refpanel,refp_df_0.8_maf_0_.001perct[,1])
dim(filtered_raw.ref_0.8_0_.001perct)
filtered_raw.seq_0.8_0_.001perct<-filter_raw(df.raw_seq,refp_df_0.8_maf_0_.001perct[,1])
dim(filtered_raw.seq_0.8_0_.001perct)

calculate_kappa(filtered_raw.ref_0.8_0_.001perct,filtered_raw.seq_0.8_0_.001perct)

####################################################################################
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
[1] dplyr_0.7.6 plyr_1.8.4

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.19      crayon_1.3.4      assertthat_0.2.0  R6_2.3.0
 [5] magrittr_1.5      pillar_1.3.0      rlang_0.2.2       data.table_1.11.8
 [9] bindrcpp_0.2.2    glue_1.3.0        purrr_0.2.5       compiler_3.4.2
[13] pkgconfig_2.0.2   bindr_0.1.1       tidyselect_0.2.5  tibble_1.4.2

"