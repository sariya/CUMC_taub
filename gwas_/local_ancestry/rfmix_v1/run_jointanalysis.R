#!/bin/Rscript

#
#Date 02/22/2019
#Sanjeev Sariya
#joint analysis per SNP adjusting for dosage and individual ancestry
#
library(data.table)
library(dplyr)
library(argparse)

#We we will use SNP dosage from imputation. 
#RFmix's ancestry component - per SNP #Run analysis per ancestry NAT, CEU and YRI.
#RFMIX_datapreparations/ancestry_components/CHR21/merged_CHR21_NAT.txt

#


parser <- ArgumentParser(description="run local ancestry for NAT and YRI adjusting for global")
parser$add_argument('-a',"--anc",help="input ancestry components per SNPs",required=TRUE)  # input components file
parser$add_argument('-d',"--dos",help="input dosage per SNPs",required=TRUE)  # input components file
args <- parser$parse_args() #make it a data structure

phenotype.file<-"input_pheno_global_ancestry.txt" # has global components too

snp_dosage.file<-normalizePath(args$dos) #viterbi file
ancestry_component.file<-normalizePath(args$anc) 

#snp_dosage.file<-"test_dosage_convert_CHR22" #this is comma separated file 
#ancestry_component.file<-"test_NAT_CHR22" #this is comma separated file

df.phenotype<-read.table(phenotype.file,header=TRUE)
df.dosage<-data.table::fread(snp_dosage.file, showProgress = TRUE,sep=",",header=FALSE)  #no header by default. Header =FALSE
df.ancestry<-data.table::fread(ancestry_component.file, showProgress = TRUE,sep=",",header=FALSE)  #no header by default. Header =FALSE

print(dim(df.dosage))
print(dim(df.ancestry))

df.dosage<-df.dosage[,-c(2:3)] #delete  a1 and a2 columns
df.ancestry<-df.ancestry[,-c(2:3)] #delete  a1 and a2 columns

if(nrow(df.dosage)!=nrow(df.ancestry)){
	stop("There is mismatch of snp counts!")
}

store_pvaluePerSNP<-as.data.frame(matrix(NA,nrow=nrow(df.dosage), ncol=3)) #per SNP P-value for V1, V3, age, sex, NAT global and YRI global 

for (i in 1:nrow(df.ancestry)){
	#
	#loop over each SNP per ancestry CEU/NAT/YRI
	df.temp_dosage<-df.dosage[i,] #get one SNP only
	snp_name<-df.temp_dosage[1,1] #store snp name 
	snp_name<-unlist(snp_name) #this is such a pain when playing data from data.table !!! :-/
	df.temp_ancestry<-df.ancestry[which(df.ancestry$V1==snp_name),] #make sure the SNP row is found

	df.merged_dosage.anc<-(as.data.frame(cbind((t((df.temp_dosage)))[-1,], (t((df.temp_ancestry)))[-1,])))

	phenosnp_cov<-cbind(df.merged_dosage.anc,df.phenotype)

	tryCatch({
		temp.df_pvalue<-summary(glm(AD  ~ as.numeric(as.character(V1)) + as.numeric(as.character(V2)) + as.numeric(as.character(age)) + as.factor(sex) +as.numeric(as.character(NAT)), 
		data = phenosnp_cov, family = "binomial"))$coeff

	},
	error=function(error_message){
		message(error_message)
		return(NA)
	})
	#catch ends for glm pvalue
	tryCatch({
	store_pvaluePerSNP[i,1]<-snp_name
	store_pvaluePerSNP[i,2]<-temp.df_pvalue[2,4] # dosage
	store_pvaluePerSNP[i,3]<-temp.df_pvalue[3,4] # local ancestry
	},
	error=function(error_message){
		message(error_message)
		return(NA)
	})
	#catch ends for storing pvalue

}#for loop ends
colnames(store_pvaluePerSNP)<-c("RSid","dosage_pvalue","ancestry_pvalue")
print(dim(store_pvaluePerSNP))

write.table(store_pvaluePerSNP, file = "", append = FALSE, quote = FALSE, sep = "\t", eol = "\n",
na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

print("Analysis ended")

