#!/bin/Rscript

#Date 02/19/2019
#Sanjeev Sariya 
#PI - Giuseppe Tosto
#use output
#
#

#Date - input - per chromsome . SNP file with names and derived allele per person per snp.
#Local ancestry is home-made from RFmix viterbi file. 
#Total columns would be N_people multipleid by 3. Three - CEU/YRI/NAT
#

library(data.table)
library(dplyr)
snp_position.file<-"rsids_allchrs_positions"
ancestry.file<-"components_CHR22_rfmix.txt" # this is space sep
snp.file<-"CHR22_snps_rfmix" #this is comma sep 
phenotype.file<-"input_pheno_global_ancestry.txt"

df.haplotype_ancestry<-data.table::fread(ancestry.file, showProgress = TRUE)  #no header by default. Header =FALSE
df.snps<-data.table::fread(snp.file,sep=",", showProgress = TRUE,header=FALSE)  #no header by default. Header =FALSE
df.phenotype<-data.table::fread(phenotype.file, showProgress = TRUE,header=TRUE)  #no header by default. Header =FALSE

df.snppos_rsids<-data.table::fread(snp_position.file, showProgress = TRUE,header=FALSE)  #no header by default. Header =FALSE

#--get CHR and positions
temp<-strsplit(df.snppos_rsids$V1,"\\:") #split column for RS 
rs_chrmatrix <- do.call(rbind,temp) #do binding 

rs_chrmatrix <- as.data.frame(rs_chrmatrix,stringsAsFactors = FALSE)
colnames(rs_chrmatrix) <- c("CHR","BP")
rs_chrmatrix$BP<-as.numeric(rs_chrmatrix$BP)
rs_chrmatrix$CHR<-as.numeric(rs_chrmatrix$CHR)
rs_chrmatrix$rsid<-df.snppos_rsids$V2

print(dim(df.haplotype_ancestry))
print(dim(df.phenotype))
print(dim(df.snps))

store_pvaluePerSNP<-as.data.frame(matrix(NA,nrow=nrow(df.snps), ncol=7)) #per SNP P-value for V1, V3, age, sex, NAT global and YRI global 
colnames(store_pvaluePerSNP)<-c("SNP_name","Pvalue_NAT.localancestry","Pvalue_YRI.localancestry","Age","Sex","NAT_globalancestry","YRI_globalancestry")
tempdf<-""
#
#Take one SNP at a time. Break it's data into per person. Then remove CEU column
#
for (i in 1:nrow(df.snps)){
print(paste(df.snps[i,1] ,i,sep="    "))
##print(df.haplotype_ancestry[i,])
print(dim(df.haplotype_ancestry[i,]))

#make df per SNP
temp_snpmatrix<-as.data.frame(matrix(df.haplotype_ancestry[i,], ncol=3, byrow=TRUE)) ###https://stackoverflow.com/questions/26973029/split-one-row-after-every-3rd-column-and-transport-those-3-columns-as-a-new-row

print(dim(temp_snpmatrix))
print(head(temp_snpmatrix))

##First column is NAT, #
#2nd is CEU,  #3rd is YRI column in haplotype SNP matrix#

#remove CEU column. We work with only two columns nows
temp_snpmatrix<-temp_snpmatrix[,-c(2)] 
print(dim(temp_snpmatrix))

phenosnp_cov<-cbind(temp_snpmatrix,df.phenotype)
print(head(phenosnp_cov))
print(class(phenosnp_cov))
tempdf<-phenosnp_cov

phenosnp_cov$unlist_V1<-unlist((phenosnp_cov$V1))
phenosnp_cov$unlist_V3<-unlist((phenosnp_cov$V3))

temp.df_pvalue<-summary(glm(AD  ~ unlist_V1 + unlist_V3 + age + sex +NAT+ YRI, data = phenosnp_cov, family = "binomial"))$coeff

store_pvaluePerSNP[i,1]<-df.snps[i,1]
store_pvaluePerSNP[i,2]<-temp.df_pvalue[2,4] # NAT local pvalue
store_pvaluePerSNP[i,3]<-temp.df_pvalue[3,4] # YRI local pvalue
store_pvaluePerSNP[i,4]<-temp.df_pvalue[4,4] #age
store_pvaluePerSNP[i,5]<-temp.df_pvalue[5,4] #sex
store_pvaluePerSNP[i,6]<-temp.df_pvalue[6,4] #NAT global pvalue
store_pvaluePerSNP[i,7]<-temp.df_pvalue[7,4] #YRI global pvalue

break
}

print(head(store_pvaluePerSNP))
pvalue_snppos_rsids<-left_join(store_pvaluePerSNP,rs_chrmatrix,by=c("SNP_name"="rsid"))

tempdf$unlist_V1<-unlist((tempdf$V1))
tempdf$unlist_V3<-unlist((tempdf$V3))

tempdf$unlist_V1<-as.factor(tempdf$unlist_V1)
tempdf$unlist_V3<-as.factor(tempdf$unlist_V3)

summary(glm(AD  ~ unlist_V1 + unlist_V3 + age + sex, data = tempdf, family = "binomial"))










