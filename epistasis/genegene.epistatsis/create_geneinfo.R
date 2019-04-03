#!/bin/usr/Rscript

#
#Date 04/03/2019
#Sanjeev Sariya
#

library("argparse")
parser <- ArgumentParser(description="Create genes info file")

parser$add <- argument('-o',"--odir",help="Store output directory",required=TRUE) #store output directory
parser$add <- argument('-x',"--pre",help="Store prefix for output file",required=TRUE) #store prefix for output file
parser$add <- argument('-b',"--bim",help="Location for bim file file",required=TRUE) ##bim file where SNPs from interested chrs are extracted
parser$add <- argument('-g',"--gene",help="Location gene file",required=TRUE) ##Comma separated gene file 

args <- parser$parse <- args() #make it a data structure

prefix.out<-args$pre ##get output prefix
file.bim <- normalizePath(args$bim) #make into full path to store 
file.gene <- normalizePath(args$gene) #make into full path and store it
out_dir<-paste(normalizePath(args$odir),"/",sep="")

##file.bim<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/EPISTASIS/ss_epistasis/genebased_otherpackages/outputagain/chromosomegenes.bim" #bim file where SNPs are present
##file.gene<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/EPISTASIS/ss_epistasis/FRGEpistasis_03122019/gene_fix_igap.csv" ##use this gene file 

print(file.bim)
print(file.gene)

df.bimfile<-read.table(file.bim,header=FALSE) 
df.genefile<-read.table(file.gene,header=TRUE,sep=",") #gene file is comma sep
 
df.bimfile<-df.bimfile[,-c(3,5,6)] ##trim it to your needed data

print(dim(df.bimfile))
print(dim(df.genefile))

number_snpswithin_gene=0
for(i in 1:nrow(df.genefile)){
    number_snpswithin_gene<-number_snpswithin_gene+length(which(df.bimfile$V4 >= df.genefile[i,3] & df.bimfile$V4 <= df.genefile[i,4] & df.bimfile$V1==df.genefile[i,2]  ))
}
##use the numbersnps with gene to create an empty data frame.

print(number_snpswithin_gene)

df.genes.info<-as.data.frame(matrix(nrow=number_snpswithin_gene,ncol=4))
iterator_snpstore=1 ##use this to iterate for snp-gene storage

##begin for loop for genes 
for(i in 1:nrow(df.genefile)){

    get_index.length<-length(which(df.bimfile$V4 >= df.genefile[i,3] & df.bimfile$V4 <= df.genefile[i,4] & df.bimfile$V1==df.genefile[i,2]  )) ##store length of SNPs found
    get_index<-(which(df.bimfile$V4 >= df.genefile[i,3] & df.bimfile$V4 <= df.genefile[i,4] & df.bimfile$V1==df.genefile[i,2]  )) ###store indices of SNPs found


    if(get_index.length>=2){
        ##If index length is more than or equal 2
        
        for(x in 1:get_index.length){
            ##Iterate over indices length. For each SNP store Gene name where it falls, it's location and Name .
            
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
        ##If length of SNPs within genomic boundary is less than 2. Print their Name
        print(df.genefile[i,])
        print(get_index.length)
    }
    ##length check ends 
    
}
##for loop ends for gene interation ends 

colnames(df.genes.info)<-c("Chromosome","Genenames","SNPnames","Position") ###set colnames 
print(dim(df.genes.info))
print(iterator_snpstore) 

df.genes.info<-df.genes.info[complete.cases(df.genes.info),] ##remove rows with no info
print(df.genes.info)


stop("we are doing good")
file.genes_info<-
write.table (df.genes.info, file = file.genes_info, append = FALSE, quote = FALSE, sep = "\t", eol = "\n",
na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
    
print("We are exiting code!!!")


