#!/bin/Rscript

#
#Date 02/21/2019
#Sanjeev Sariya 
#Take in coded RFmix output and calculate per person switch

#
#Column are persons and rows are SNPs #take in Viterbi output
library(data.table)
library(argparse)
parser <- ArgumentParser(description="calculate switches using haps ancestry rfmix v1 viterbi output")
parser$add_argument('-v',"--vit",help="Number",required=TRUE)
parser$add_argument('-x',"--pre",help="prefix for output file",required=TRUE)

args <- parser$parse_args() #make it a data structure

rfmix_coded.file<-normalizePath(args$vit) #viterbi file

df.rfmixcoded<-data.table::fread(rfmix_coded.file, showProgress = TRUE,header=FALSE)
print(dim(df.rfmixcoded))

num_snps=nrow(df.rfmixcoded)
num_people=ncol(df.rfmixcoded)/2

print(paste("the number of SNPs are",num_snps))
print(paste("the number of people are",num_people))

df_countperperson <- data.frame(counts=matrix(integer(), ncol = 1, nrow = num_people),stringsAsFactors=FALSE) #store switch counts in it

itr_person=0 #use it as iterator

##loop per person first then loop within that person using SNPs data#
while(itr_person<num_people){
	itr_person<-itr_person+1

	range_col=itr_person*2 # upto this person's column

	range_col<-as.integer(as.numeric(range_col))

	df_perperson<-df.rfmixcoded[,c(range_col-1,range_col),with=FALSE] #https://stackoverflow.com/questions/34798957/subset-by-column-index-in-r-data-table-vs-dataframe
	count_switch_perperson=0	#use this per person and get counts
	
	#iterate over SNPs now
	for(i in 1:(num_snps-1)){

		if(df_perperson[i,1]!=df_perperson[i+1,1] | df_perperson[i,2]!=df_perperson[i+1,2]){
			count_switch_perperson=count_switch_perperson+1
		}
	}
	print(paste("Done with person",itr_person))
	#for loop ends per person

	df_countperperson[itr_person,1]<-count_switch_perperson #set number of switch in null dataf
}
#while loop ends

print(dim(df_countperperson))

write.table(x, file =args$pre , append = FALSE, quote = FALSE, sep = " ", eol = "\n",
na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape","double"), fileEncoding = "")

print("Exiting code")


