#!/usr/bin/Rscript

#
#Date 05/10/2019
#Sanjeev Sariya
#05/10/2019

#
#Input file is viterbi output from RFmix. It has twice the number of individuals column counts
# 
library("argparse")

#
#RFmix 1.5.X
#The output has admixed individuals columns. Per person one colum
#1- is NAT, 2 is CEU and 3- is Yurobian as the information in output  #--it depends how your input in class file is designed.  

##

parser <- ArgumentParser(description="make ancetry from haps in rfmix v1 output")
parser$add_argument('-v',"--vit",help="input viterbifile",required=TRUE) #
parser$add_argument('-x',"--outpre",help="output prefix",required=TRUE) # prefix for output
args <- parser$parse_args() #make it a data structure

viterbi.file<-args$vit
prefix<-args$outpre

##file.input<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/output_inputfiles_RFmix/CHR22/CHR22_rfmix.0.Viterbi.txt"
##df.viterbi<-read.table(file.input,header=FALSE)
print(viterbi.file)
print(prefix)

df.viterbi<-read.table(viterbi.file,header=FALSE)
print(dim(df.viterbi))

##use lapply .. khi-khi

print(nrow(df.viterbi))
print(ncol(df.viterbi))

#
#############################################
#get_ancestry function starts
#############################################
get_ancestry<-function(x,y){
#1- is NAT, 2 is CEU and 3- is Yurobian

ancestry_mat<-matrix(0,nrow = 1, ncol = 1)

#--if 100% from one ancestry
if(x==y){
    if(x==1){
        ancestry_mat[1,1]<-1
}
#1- is NAT, 2 is CEU and 3- is Yurobian
if(x==2){
    ancestry_mat[1,1]<-2
}
#1- is NAT, 2 is CEU and 3- is Yurobian
if(x==3){
    ancestry_mat[1,1]<-3
}
#1- is NAT, 2 is CEU and 3- is Yurobian
} else{ 

#--if 100% not from one ancestry

if((x==1 & y==2) |(x==2 & y==1) ){
    ancestry_mat[1,1]<-4
}
if ((y==2 & x==3) | (x==2 & y==3)){
    ancestry_mat[1,1]<-5
}

#1- is NAT, 2 is CEU and 3- is Yurobian
if((x==1 & y==3) | ( x==3 & y==1 )){
    ancestry_mat[1,1]<-6
}

} ##else ends 

return(ancestry_mat)
}
#############################################
#get_ancestry function ends
#############################################


haploStates_ancestry_per_person<-apply(df.viterbi,1,function(i)
{

#--make data frame of row
row_df<-t(as.data.frame(i))

iterate<-ncol(row_df)/2 #get variable to iterate over columns in a row

#--iterate over columns
ancestral_components<-sapply(1:iterate,function(index){
	
    t<-index*2 #t is used for index in a vector
	
    output_mat<-get_ancestry(row_df[t-1],row_df[t])
    return(output_mat)
})
      
})

print(dim(haploStates_ancestry_per_person))
print(dim(t(haploStates_ancestry_per_person)))

outfile<-paste(prefix,".txt",sep="")
write.table(t(haploStates_ancestry_per_person), file = outfile, append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE) #write to an output file

print("check current working directory")



