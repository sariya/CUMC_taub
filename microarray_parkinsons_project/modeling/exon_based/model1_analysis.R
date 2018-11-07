#!/bin/Rscript

#--Date: 06/04/2018
#Sanjeev Sariya Dr. Tosto 
#--Karen Mardar PD expression data
library("lme4",lib.loc="/home/ss5505/R_LIB/")

library(dplyr)
library("argparse")
parser <- ArgumentParser(description="PD expression analysis")

parser$add_argument('-e',"--exon",help="Location for exon chunk file",required=TRUE)
parser$add_argument('-n',"--num",help="Number",required=TRUE)
parser$add_argument('-o',"--out",help="output directory",required=TRUE)
parser$add_argument('-x',"--pre",help="prefix for output file",required=TRUE)

args <- parser$parse_args() #make it a data structure
print(normalizePath(args$out))
out<-normalizePath(args$out)
chunk<-args$num
exon_file<-args$exon
print(exon_file)
prefix<-args$pre
exon.dataf<-read.table(exon_file,sep="\t",header=TRUE)
print(dim(exon.dataf))
#-------

row.names(exon.dataf)<-exon.dataf$id #--add row names
transpose_data.exon<-t(exon.dataf)

pheno<-read.table("pheno_rin_study",sep="\t",header=TRUE) # Pheno with Study and RIN
print(head(pheno))
row.names(pheno)<-pheno$IID
PDmerge<-merge(transpose_data.exon,pheno, by=0,all=TRUE)
print("merge")
PDmerge<-PDmerge[-c(1),] #--delete useless rows
PDmerge<-PDmerge[,-c(1)] #delete x column
names<-colnames(PDmerge)
print("names done")

PDmerge$AGE_SCALE<-(scale(PDmerge$AGE))[,1]
missing_celpheno<-which(is.na(PDmerge$IID))
PDmerge<-PDmerge[-c(missing_celpheno),]

exon_mod<-as.data.frame(matrix(NA,nrow=ncol(PDmerge)-10, ncol=2))
#--store pvalue and corresponding transcript id

print("Begin loop modeling")

for(i in 1:(ncol(PDmerge)-10) ){

tryCatch({
	exon_mod[i,1] <- as.numeric(as.character(names[i])) #  storeexon
},
	error=function(error_message){
	message(error_message)
	return(NA)
}
	)
	tryCatch(
	{

	exon_mod[i,2]<-summary( glmer(PD~as.factor(SEX)+AGE_SCALE+as.numeric(as.character(PDmerge[,i]))+(1|FID) + (1|STUDY) + RIN , data=PDmerge,family = binomial(link=logit)) )$coeff[4,4]
	},
	error=function(error_message){
	message(error_message)
	return(NA)
	}
	)
}
#--
print("Done with Loop modeling")

colnames(exon_mod)<-c("Probeset","Pvalue")

temp_file<-paste(args$pre,"_output",sep="")
out_file<-paste(out,temp_file,sep="/")
print(out_file)
write.table(exon_mod,out_file,quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
print("Check output file please")