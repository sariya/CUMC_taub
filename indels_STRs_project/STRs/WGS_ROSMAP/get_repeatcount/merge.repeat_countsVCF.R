#!/bin/Rscript

#
#Date 05/22/2019
#Sanjeev Sariya
##badri, Giuseppe 

##/usr/local/bin/Rscript
#
#Read in VCF#Create genotype based on repeat count and allele information 
#from VCF GT : 0|0, 0|1, so on and so forth#

#
#Rscript  merge.repeat_countsVCF.R  -v VCF -r REP -o ODIR -p PREFIX
#
library(vcfR) ##This is vcfR 1.8.0
library(argparse)

parser <- ArgumentParser(description="create genotype from VCF and motif repeat files") #FAST and GCTA --
parser$add_argument('-v',"--vcf",help="Location for vcf file",required=TRUE) # 
parser$add_argument('-r',"--rep",help="Location for vcf file",required=TRUE) #
parser$add_argument('-o',"--odir",help="Store output directory",required=TRUE) #store output directory
parser$add_argument('-p',"--prefix",help="Store output prefix",required=TRUE) #store output prefix

args <- parser$parse_args() #make it a data structure
prefix<-args$prefix
out_dir<-normalizePath(args$odir)
out.file<-paste(out_dir,prefix,sep="/")

##file.repeat_count<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/calculate_dosage.byrepeats/repeatcounts_CHR18"

file.repeat_count<-normalizePath(args$rep)
df.repeat_count<-read.table(file.repeat_count,header=FALSE)
print(dim(df.repeat_count))

#vcf_file<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/cleaning_mergedSTRs/cleanedmiss_depth_extractedpersons.final/filter.str.no_monomorphic.poor_depth.highmissing.postcleaning.removeindi_CHR18.recode.vcf"

vcf_file<-normalizePath(args$vcf)
vcf <- read.vcfR( vcf_file, verbose = TRUE )

name_strs<-as.data.frame(matrix(nrow=nrow(vcf@gt),ncol=4))
name_strs[,1]<-vcf@fix[,3] ##set STR names
name_strs[,2]<-paste(vcf@fix[,1] ,vcf@fix[,2],vcf@fix[,3],sep=":")  ##get STR names

person.names<-names(vcf@gt[1,2:ncol(vcf@gt)])
store.genotype<-(matrix(nrow=nrow(vcf@gt),ncol=ncol(vcf@gt)-1)) ##per person

###for(i_snpnumber in 1:2 ){ 

for(i_snpnumber in 1:nrow(vcf@gt) ){ 

    index_str=which(df.repeat_count$V1==name_strs[i_snpnumber,1])
    
    if(length( index_str)!=0){
        
        temp.str.info<-vcf@gt[i_snpnumber,] 

        for(x_personnum in 2:length(names(temp.str.info))){
            ##loop per person      ##1st position value is "FORMAT"

            get.genotype<- unlist(strsplit(temp.str.info[x_personnum], ":", fixed = FALSE, perl = FALSE, useBytes = FALSE))[1] 

            if(grepl("|", get.genotype,fixed=TRUE)==TRUE){
                
                getref_alt<-unlist(strsplit(get.genotype,"|", fixed = FALSE, perl = FALSE, useBytes = FALSE))
                
                if(getref_alt[1]!="0" & getref_alt[3]!="0"){
                    
                    repeats_alteranteallele<-df.repeat_count[index_str,3] ##store count value
                    
                    ##print(repeats_alteranteallele)					print("multiple alleles in alterante ")
                    ##print(class(repeats_alteranteallele)) 					print(unlist(strsplit( as.character(repeats_alteranteallele) , ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)) )			
                    ##split and unlist values. If comma present good, if not. no worries
                    
                    counts.allelemotifs<-(unlist(strsplit( as.character(repeats_alteranteallele) , ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)) ) 

                    ##print(length(counts.allelemotifs)) 					print(counts.allelemotifs[as.numeric(getref_alt[1])])
                    ##print(counts.allelemotifs[as.numeric(getref_alt[3])]) 					print(as.numeric(counts.allelemotifs[as.numeric(getref_alt[1])])+ as.numeric(counts.allelemotifs[as.numeric(getref_alt[3])]))
                    store.genotype[i_snpnumber,x_personnum-1]<- as.numeric(counts.allelemotifs[as.numeric(getref_alt[1])])+ as.numeric(counts.allelemotifs[as.numeric(getref_alt[3])]) 
                    
                } else if(getref_alt[1]=="0" & getref_alt[3]=="0"){
                    
                    store.genotype[i_snpnumber,x_personnum-1]<-df.repeat_count[index_str,2]*2
                }else if(getref_alt[1]!="0" & getref_alt[3]=="0"){
                    
                    repeats_alteranteallele<-df.repeat_count[index_str,3] ##store count value
                    counts.allelemotifs<-(unlist(strsplit( as.character(repeats_alteranteallele) , ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)) )
                    store.genotype[i_snpnumber,x_personnum-1]<- as.numeric(counts.allelemotifs[as.numeric(getref_alt[1])]) + df.repeat_count[index_str,2]
                    
                } else if(getref_alt[1]=="0" & getref_alt[3]!="0"){
                    
                    repeats_alteranteallele<-df.repeat_count[index_str,3] ##store count value
                    counts.allelemotifs<-(unlist(strsplit( as.character(repeats_alteranteallele) , ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)) )
                    store.genotype[i_snpnumber,x_personnum-1]<- as.numeric(counts.allelemotifs[as.numeric(getref_alt[3])]) + df.repeat_count[index_str,2]
                }
            }
            ##if | found in genotype
            else{
                
                if(get.genotype=="."){
                    
                    store.genotype[i_snpnumber,x_personnum-1]<-NA ##fix here for xperson
                    
                } else if(grepl("/", get.genotype,fixed=TRUE)==TRUE & get.genotype=="./."){
                    
                    store.genotype[i_snpnumber,x_personnum-1]<-NA ##fix here for xperson
                } else{
                    
                    print(get.genotype)
                    stop("We do not know genotype for this")
                }                
                
            }
            ##no | in genotype value ends
        }
        ##for loop ends for per person genotype
    }
    ##if length ==0 ends
    else{
        print(paste(name_strs[i_snpnumber,1]," missing from counted motifs"))
    }
    ##if else check ends that STR was found
}
##for loop ends

colnames(store.genotype)<-person.names
colnames(name_strs)<-c("STR_name","harmonized_STR_name","allele1","allele2")
name_strs[,3:4]<-"AABBCC"
strnames.genotype<-cbind(name_strs,store.genotype)
strnames.genotype<-strnames.genotype[,-c(1)]

print(paste(out.file,"withoutcolnames",sep="_"))
print(paste(out.file,"withcolnames",sep="_"))

write.table (strnames.genotype, file = paste(out.file,"withoutcolnames",sep="_"), append = FALSE, quote = FALSE, sep = ",", eol = "\n",
             na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape","double"), fileEncoding = "")

write.table (strnames.genotype, file = paste(out.file,"withcolnames",sep="_"), append = FALSE, quote = FALSE, sep = ",", eol = "\n",
             na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape","double"), fileEncoding = "")

print("Exiting code")
