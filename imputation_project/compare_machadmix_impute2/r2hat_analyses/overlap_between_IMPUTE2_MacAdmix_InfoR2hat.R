#!/bin/Rscript

#
#Find overlap between two SNP list
#
#Get SNPs in different info and rsq threshold and MAF bins
# Calculate mean for those snps
#Imputation from 1000G
#Sanjeev Sariya
#10/24/2018
#

# cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/accuracy_imputation/r2_hat_calculation/machadmix/calculate_r2hat_perchunk/10242018_calculation

library(dplyr)
file.merged_machadmix<-"joined_r2hat_machadmixinfo_CHR21_cleaned"
file.merged_impute2<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/accuracy_imputation/r2_hat_calculation/impute2/1000G/calculate_r2hat_perchunk/final_probs/joined_r2hat_imputeinfo_CHR21_cleaned"

df.machadmix_merged<-data.table::fread(file.merged_machadmix, showProgress = TRUE,header=TRUE)
df.impute2_merged<-data.table::fread(file.merged_impute2, showProgress = TRUE,header=TRUE)

print(dim(df.machadmix_merged))
print(dim(df.impute2_merged))

#https://stackoverflow.com/questions/3695677/how-to-find-common-elements-from-multiple-vectors

print(length(Reduce(intersect, list(df.impute2_merged$rs_id,df.machadmix_merged$SNP))))
#549091

print(length(Reduce(intersect, list(df.machadmix_merged$SNP,df.impute2_merged$rs_id))))
#549091
overlappingSNPs<-Reduce(intersect, list(df.machadmix_merged$SNP,df.impute2_merged$rs_id))

#get dimensions
print(dim(as.data.frame(overlappingSNPs)))

#write to a outfile with a header
write.table(as.data.frame(overlappingSNPs),"overlappingSNPs_header",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

#make data frame
overlappingSNPs<-as.data.frame(overlappingSNPs)

#get Mach SNps and merged info and R2hat
machadmix_joined_overlapping_r2Info<-left_join(overlappingSNPs,df.machadmix_merged,by=c("overlappingSNPs"="SNP"))

#get IMPUTE2 SNps and merged info and R2hat
impute2_joined_overlapping_r2Info<-left_join(overlappingSNPs,df.impute2_merged,by=c("overlappingSNPs"="rs_id"))

write.table(impute2_joined_overlapping_r2Info,"overlap_IMPUTE2InfoR2Hat",quote=FALSE, col.names=TRUE, row.names=FALSE,sep="\t")
write.table(machadmix_joined_overlapping_r2Info,"overlap_machadmixInfoR2Hat",quote=FALSE, col.names=TRUE, row.names=FALSE,sep="\t")

mean(machadmix_joined_overlapping_r2Info$R2_hat)
mean(impute2_joined_overlapping_r2Info$R2_hat)


#
#Info >=0.40
#
#compute for IMPUTE2 overllaping SNPs
impute2_joined_overlapping_r2Info_0.40<-impute2_joined_overlapping_r2Info[impute2_joined_overlapping_r2Info$info>=0.40,]
dim(impute2_joined_overlapping_r2Info_0.40) ## 371108     16

min(impute2_joined_overlapping_r2Info$info)
max(impute2_joined_overlapping_r2Info$info)
max(impute2_joined_overlapping_r2Info$exp_freq_a1)
min(impute2_joined_overlapping_r2Info$exp_freq_a1)

impute2_joined_overlapping_r2Info_0.40_1_5MAF<-impute2_joined_overlapping_r2Info_0.40[which(impute2_joined_overlapping_r2Info_0.40$exp_freq_a1 >=0.01 & impute2_joined_overlapping_r2Info_0.40$exp_freq_a1 <=0.05) ,]
print( dim(impute2_joined_overlapping_r2Info_0.40_1_5MAF)) #79679
print(mean(impute2_joined_overlapping_r2Info_0.40_1_5MAF$R2_hat)) # 0.92

impute2_joined_overlapping_r2Info_0.40_0.1_1MAF<-impute2_joined_overlapping_r2Info_0.40[which(impute2_joined_overlapping_r2Info_0.40$exp_freq_a1 >=0.001 & impute2_joined_overlapping_r2Info_0.40$exp_freq_a1 <0.01) ,]
print( dim(impute2_joined_overlapping_r2Info_0.40_0.1_1MAF)) # 191864
print(mean(impute2_joined_overlapping_r2Info_0.40_0.1_1MAF$R2_hat)) # 0.8309268


impute2_joined_overlapping_r2Info_0.40_0_0.1MAF<-impute2_joined_overlapping_r2Info_0.40[which(impute2_joined_overlapping_r2Info_0.40$exp_freq_a1 >=0& impute2_joined_overlapping_r2Info_0.40$exp_freq_a1 <0.001 ) ,]
print( dim(impute2_joined_overlapping_r2Info_0.40_0_0.1MAF)) # 9963
print(mean(impute2_joined_overlapping_r2Info_0.40_0_0.1MAF$R2_hat)) # 0.6250607

#
#Rsq >=0.30
#
#compute for Mach dmix overlapping SNPs

machadmix_joined_overlapping_r2Info_0.30<-machadmix_joined_overlapping_r2Info[machadmix_joined_overlapping_r2Info$Rsq >=0.30,]
print(dim(machadmix_joined_overlapping_r2Info_0.30))
print(min(machadmix_joined_overlapping_r2Info_0.30$Rsq))
print(max(machadmix_joined_overlapping_r2Info_0.30$Rsq))

machadmix_joined_overlapping_r2Info_0.30_1_5MAF<-machadmix_joined_overlapping_r2Info_0.30[which(machadmix_joined_overlapping_r2Info_0.30$MAF.x >=0.01& machadmix_joined_overlapping_r2Info_0.30$MAF.x<=0.05),]
print(dim(machadmix_joined_overlapping_r2Info_0.30_1_5MAF))
print(mean(machadmix_joined_overlapping_r2Info_0.30_1_5MAF$R2_hat))

machadmix_joined_overlapping_r2Info_0.30_0.1_1MAF<-machadmix_joined_overlapping_r2Info_0.30[which(machadmix_joined_overlapping_r2Info_0.30$MAF.x >=0.001 & machadmix_joined_overlapping_r2Info_0.30$MAF.x<0.01),]
print(dim(machadmix_joined_overlapping_r2Info_0.30_0.1_1MAF))
print(mean(machadmix_joined_overlapping_r2Info_0.30_0.1_1MAF$R2_hat))

machadmix_joined_overlapping_r2Info_0.30_0_0.1MAF<-machadmix_joined_overlapping_r2Info_0.30[which(machadmix_joined_overlapping_r2Info_0.30$MAF.x >=0& machadmix_joined_overlapping_r2Info_0.30$MAF.x <0.001),]
print(dim(machadmix_joined_overlapping_r2Info_0.30_0_0.1MAF))
print(mean(machadmix_joined_overlapping_r2Info_0.30_0_0.1MAF$R2_hat))




