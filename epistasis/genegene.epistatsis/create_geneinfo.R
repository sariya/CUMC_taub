#!/bin/usr/Rscript

#
#Date 04/03/2019
#Sanjeev Sariya
#


file.bim<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/EPISTASIS/ss_epistasis/genebased_otherpackages/outputagain/chromosomegenes.bim"
file.gene<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/EPISTASIS/ss_epistasis/FRGEpistasis_03122019/gene_fix_igap.csv"

df.bimfile<-read.table(file.bim,header=FALSE)
df.genefile<-read.table(file.gene,header=TRUE,sep=",")
 
df.bimfile<-df.bimfile[,-c(3,5,6)] ##trim it to your needed data

print(dim(df.bimfile))
print(dim(df.genefile))


number_snpswithin_gene=0
for(i in 1:nrow(df.genefile)){
number_snpswithin_gene<-number_snpswithin_gene+length(which(df.bimfile$V4 >= df.genefile[i,3] & df.bimfile$V4 <= df.genefile[i,4] & df.bimfile$V1==df.genefile[i,2]  ))
}

print(number_snpswithin_gene)

df.genes.info<-as.data.frame(matrix(nrow=number_snpswithin_gene,ncol=4))
iterator_snpstore=1 ##use this to iterate for snp-gene storage

##begin for loop for genes 
for(i in 1:nrow(df.genefile)){

get_index.length<-length(which(df.bimfile$V4 >= df.genefile[i,3] & df.bimfile$V4 <= df.genefile[i,4] & df.bimfile$V1==df.genefile[i,2]  ))
get_index<-(which(df.bimfile$V4 >= df.genefile[i,3] & df.bimfile$V4 <= df.genefile[i,4] & df.bimfile$V1==df.genefile[i,2]  ))


if(get_index.length>=2){

for(x in 1:get_index.length){

temp.snpname<-as.character(df.bimfile[get_index[x],2])
temp.snpposition<-df.bimfile[get_index[x],3]
temp.genename<-as.character(df.genefile[i,1])

df.genes.info[iterator_snpstore,1]<-df.genefile[i,2] #store chromosome
df.genes.info[iterator_snpstore,2]<-temp.genename #store gene name
df.genes.info[iterator_snpstore,3]<-temp.snpname #store SNPnames 
df.genes.info[iterator_snpstore,4]<-temp.snpposition #store SNP position 

iterator_snpstore<-iterator_snpstore+1 ##increase the iterator
}

##for loop for index loop ends
} else{
print(df.genefile[i,])
print(get_index.length)
}
##length check ends 

}

##for loop ends for gene interation ends 


colnames(df.genes.info)<-c("Chromosome","Genenames","SNPnames","Position")
print(dim(df.genes.info))
print(iterator_snpstore)

df.genes.info<-df.genes.info[complete.cases(df.genes.info),]


