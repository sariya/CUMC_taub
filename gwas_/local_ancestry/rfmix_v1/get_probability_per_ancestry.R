#!/bin/Rscript


#
#Sanjeev Sariya
#PI Dr. Giuseppe Tosto
#Date 15 Dec 2018
#PH 19
#
#RFmix 1.5.X
#The output has2Xinput admixed individuals columns
#1- is NAT, 2 is CEU and 3- is Yurobian as the information in output  #--it depends how your input in class file is designed.  


library(data.table)
library("argparse")
parser <- ArgumentParser(description="make ancetry from haps in rfmix v1 output")
parser$add_argument('-v',"--vit",help="input viterbifile",required=TRUE) #
parser$add_argument('-x',"--outpre",help="output prefix",required=TRUE) # prefix for output
args <- parser$parse_args() #make it a data structure

viterbi.file<-args$vit
prefix<-args$outpre


#viterbi.file<-"againoutput.0.Viterbi.txt"
df.viterbi<-fread(viterbi.file)
print(nrow(df.viterbi))
print(ncol(df.viterbi))

#############################################
#get_ancestry function starts
#############################################
get_ancestry<-function(x,y){
#1- is NAT, 2 is CEU and 3- is Yurobian

ancestry_mat<-matrix(0,nrow = 1, ncol = 3)

#--if 100% from one ancestry
if(x==y){
if(x==1){
ancestry_mat[1,1]<-2
}
#1- is NAT, 2 is CEU and 3- is Yurobian
if(x==2){
ancestry_mat[1,2]<-2
}
#1- is NAT, 2 is CEU and 3- is Yurobian
if(x==3){
ancestry_mat[1,3]<-2
}
#1- is NAT, 2 is CEU and 3- is Yurobian
} else{ 

#--if 100% not from one ancestry

if(x==1){
ancestry_mat[1,1]<-1
}
if (y==1){
ancestry_mat[1,1]<-1
}

#1- is NAT, 2 is CEU and 3- is Yurobian
if(x==2){
ancestry_mat[1,2]<-1
}
if (y==2){
ancestry_mat[1,2]<-1
}
#1- is NAT, 2 is CEU and 3- is Yurobian
if(x==3){
ancestry_mat[1,3]<-1
}
if (y==3){
ancestry_mat[1,3]<-1
} }

return(ancestry_mat)
}
#############################################
#get_ancestry function ends
#############################################

snp_ancestry_per_person<-apply(df.viterbi,1,function(i)
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

return(ancestral_components)
}
)


print(dim(snp_ancestry_per_person))
print(dim(t(snp_ancestry_per_person)))

outfile<-paste(prefix,".txt",sep="")
write.table(t(snp_ancestry_per_person), file = outfile, append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE) #write to an output file

print("check current working directory")



