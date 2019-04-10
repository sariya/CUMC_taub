#!/bin/Rscript

#
#Date 04/09/2019
#Sanjeev Sariya
##badri, Giuseppe 

##/usr/local/bin/Rscript
#
#Read in VCF
#

library(vcfR)
library(argparse)

parser <- ArgumentParser(description="create genotype from VCF files") #FAST and GCTA --
parser$add_argument('-v',"--vcf",help="Location for vcf file",required=TRUE) # 
parser$add_argument('-o',"--odir",help="Store output directory",required=TRUE) #store output directory
parser$add_argument('-p',"--prefix",help="Store output prefix",required=TRUE) #store output prefix

args <- parser$parse_args() #make it a data structure
prefix<-args$prefix
vcf_file<-normalizePath(args$vcf)

out_dir<-normalizePath(args$odir)
out.file<-paste(out_dir,prefix,sep="/")
print(out.file)
print(vcf_file)

##vcf_file<-"filter.str.no_monomorphic.poor_depth.highmissing.postcleaning.removeindi_CHR22.recode.vcf"
vcf <- read.vcfR( vcf_file, verbose = TRUE )

name_strs<-as.data.frame(matrix(nrow=nrow(vcf@gt),ncol=3))
name_strs[,1]<-paste(vcf@fix[,1] ,vcf@fix[,2],vcf@fix[,3],sep=":")

person.names<-names(vcf@gt[1,2:ncol(vcf@gt)])
store.genotype<-(matrix(nrow=nrow(vcf@gt),ncol=ncol(vcf@gt)-1)) ##per person

for(i_snpnumber in 1:nrow(vcf@gt)){ 

###for(i_snpnumber in 1:2){

    print(paste("we'll begin",name_strs[i_snpnumber,1],sep=" ")) ####fix i here

    temp.str.info<-vcf@gt[i_snpnumber,] 

    for(x_personnum in 2:length(names(temp.str.info))){
        ##loop over each person for genotype
        
        get.genotype<- unlist(strsplit(temp.str.info[x_personnum], ":", fixed = FALSE, perl = FALSE, useBytes = FALSE))[1] 

        if(grepl("|", get.genotype,fixed=TRUE)==TRUE){
            ##if | found in genotype
            getref_alt<-unlist(strsplit(get.genotype,"|", fixed = FALSE, perl = FALSE, useBytes = FALSE))

            if(getref_alt[1]!="0" & getref_alt[3]!="0"){
                ##both non reference
                store.genotype[i_snpnumber,x_personnum-1]<-2 
            }else if(getref_alt[1]=="0" & getref_alt[3]!="0"){
                ##one reference, one non ref
                store.genotype[i_snpnumber,x_personnum-1]<-1 
            } else if(getref_alt[1]=="0" & getref_alt[3]=="0"){
                ##both reference
                store.genotype[i_snpnumber,x_personnum-1]<-0 ##fix here for xperson

            } else if(getref_alt[1]!="0" & getref_alt[3]=="0"){
                ##again one ref and one non refere
                store.genotype[i_snpnumber,x_personnum-1]<-1 
            }
       ##grepl ends for "|"
        } else{
            if(get.genotype=="."){
                ##if dot. We make it NA
                store.genotype[i_snpnumber,x_personnum-1]<-NA ##fix here for xperson

            } else if(grepl("/", get.genotype,fixed=TRUE)==TRUE & get.genotype=="./."){
                ##again make it NA
                store.genotype[i_snpnumber,x_personnum-1]<-NA ##fix here for xperson

            }else{
                ##anything else
                print(get.genotype)
                stop("We do not know genotype for this")

            } ##else ends for . match

        } ##else ends for no | 

    } ##loop ends for per person

    print(paste("Done through",name_strs[i_snpnumber,1],sep=" ")) ##fix i here
} ##iterate over all STRs

print(dim(store.genotype))

colnames(store.genotype)<-person.names
colnames(name_strs)<-c("STR_name","allele1","allele2")
name_strs[,2:3]<-"AABBCC"

strnames.genotype<-cbind(name_strs,store.genotype)

print(paste(out.file,"withoutcolnames",sep="_"))
print(paste(out.file,"withcolnames",sep="_"))

write.table (strnames.genotype, file = paste(out.file,"withoutcolnames",sep="_"), append = FALSE, quote = FALSE, sep = ",", eol = "\n",
    na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape","double"), fileEncoding = "")
	
write.table (strnames.genotype, file = paste(out.file,"withcolnames",sep="_"), append = FALSE, quote = FALSE, sep = ",", eol = "\n",
    na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape","double"), fileEncoding = "")

print("Exiting code!!")




