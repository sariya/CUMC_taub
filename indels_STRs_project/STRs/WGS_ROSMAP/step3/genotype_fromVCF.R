#!/bin/Rscript

#
#Date 04/09/2019
#Sanjeev Sariya
##badri, Giuseppe 

#
#Read in VCF
#

library(vcfR)
vcf_file<-"filter.str.no_monomorphic.poor_depth.highmissing.postcleaning.removeindi_CHR22.recode.vcf"
vcf <- read.vcfR( vcf_file, verbose = TRUE )

name_strs<-as.data.frame(matrix(nrow=nrow(vcf@gt),ncol=3))
name_strs[,1]<-paste(vcf@fix[,1] ,vcf@fix[,2],vcf@fix[,3],sep=":")

person.names<-names(vcf@gt[1,2:ncol(vcf@gt)])
store.genotype<-(matrix(0,nrow=nrow(vcf@gt),ncol=ncol(vcf@gt)-1)) ##per person

for(i_snpnumber in 1:nrow(vcf@gt)){ 

##for(i in 1:2){

print(paste("we'll begin",name_strs[i_snpnumber,1],sep=" ")) ####fix i here

temp.str.info<-vcf@gt[i_snpnumber,] 

for(x_personnum in 2:length(names(temp.str.info))){

get.genotype<- unlist(strsplit(temp.str.info[x_personnum], ":", fixed = FALSE, perl = FALSE, useBytes = FALSE))[1] ##fix here for xperson

if(grepl("|", get.genotype,fixed=TRUE)==TRUE){
getref_alt<-unlist(strsplit(get.genotype,"|", fixed = FALSE, perl = FALSE, useBytes = FALSE))

if(getref_alt[1]!="0"){

store.genotype[i_snpnumber,x_personnum-1]<-2 ##fix here for xperson

}else if(getref_alt[1]=="0" & getref_alt[3]!="0"){

store.genotype[i_snpnumber,x_personnum-1]<-1 ##fix here for xperson

} else if(getref_alt[1]=="0" & getref_alt[3]=="0"){ 

store.genotype[i_snpnumber,x_personnum-1]<-0 ##fix here for xperson

}
##grepl ends for "|"
} else{
if(get.genotype=="."){

store.genotype[i_snpnumber,x_personnum-1]<-0 ##fix here for xperson
} else if(grepl("/", get.genotype,fixed=TRUE)==TRUE){
store.genotype[i_snpnumber,x_personnum-1]<-0 ##fix here for xperson

}else{
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

write.table (strnames.genotype, file = "withoutcolnames", append = FALSE, quote = FALSE, sep = ",", eol = "\n",
    na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape","double"), fileEncoding = "")
	
write.table (strnames.genotype, file = "withcolnames", append = FALSE, quote = FALSE, sep = ",", eol = "\n",
    na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape","double"), fileEncoding = "")

print("Exiting code!!")




