#!/bin/Rscript


library(dplyr)

#
#Date 10/25/2018
#Sanjeev Sariya

#
#You need pheno file with age, sex
#You need MDS file from King
#You need order file - columbia, MESA, and etc studies

#
file.pheno<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/GWAS_data_analyses/PHENO/WHICAPEFIGA_1066PR_MES_NOMAS_finealpheno_06072018"
df.pheno<-read.table(file.pheno,header=TRUE)

file.order_study<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/GWAS_data_analyses/PHENO/MESA_HGWAS1to6PR_NOMAS_merged.sample"
df.order_study<-read.table(file.order_study,header=TRUE)
print(dim(df.order_study))

key<-read.table("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/keys/old_ids_new_ids_hgwas_pr_ordered_dosage", header=T)
print(dim(key))

#-- read sample ids from merged GWAS batches
sample_ids_merged<-read.table("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/GWAS_data_analyses/merged_gen/MESA_HGWAS1to6PR_NOMAS_merged_CHR6.sample", header=T)
head(sample_ids_merged)

#--delete row 1
sample_ids_merged<-sample_ids_merged[-c(1),]

#--only keep 1 and 2nd column
sample_ids_merged<-sample_ids_merged[,c(1,2)]

#--name them so we perform left_join later on
colnames(sample_ids_merged)<-c("OldFID","OldIID")
head(sample_ids_merged)
print(dim(sample_ids_merged))

mds_removed_indiv<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/GWAS_data_analyses/KINSHIP_06182018/gwas/gwas_nomasmesa__PRhgwas126/removed_indiv_MDS_King_output"
df.mds<-read.table(mds_removed_indiv,header=TRUE)

print(dim(df.mds))

#--join order, old ids and new ids
order_oldnew.keys<-left_join(key,df.order_study,by=c("OldIID"="IID"))
print(dim(order_oldnew.keys))

#--check where order/study column are fudged
which(is.na(order_oldnew.keys$STUDY))

#--join order, old keys to phenotype
pheno_order_oldnew.keys<-left_join(order_oldnew.keys,df.pheno,by=c("New_IID"="IID"))
print(dim(pheno_order_oldnew.keys))

mds_pheno.ordered.oldnew.keys<-left_join(pheno_order_oldnew.keys,df.mds,by=c("New_IID"="V2"))
print(dim(mds_pheno.ordered.oldnew.keys))

ordered_merged_mds_ad_order_new.oldIDS <- mds_pheno.ordered.oldnew.keys[order(mds_pheno.ordered.oldnew.keys$ORDER),] 
write.table(ordered_merged_mds_ad_order_new.oldIDS,"mergedmds_oldnewIIDs_order_10252018",row.names=FALSE, col.names=TRUE,quote=FALSE, sep="\t")


