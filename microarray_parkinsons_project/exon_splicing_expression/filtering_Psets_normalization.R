#!/bin/Rscript

#
#
#Sanjeev Sariya 11/07/2018
#Sanjeev Sariya
#
#GiusTosto KarMarder PD HuExST1.0

library("oligo",lib.loc="~/R_LIB")
library("pd.huex.1.0.st.v2",lib.loc="~/R_LIB") 

setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/CEL/")

##############################################################################
#Reading Cel files
##############################################################################
celFiles <- list.celfiles("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/CEL")
rawData <- read.celfiles(celFiles)
callps <- paCalls(rawData,method="PSDABG") 

##############################################################################
#Filtering Probesets Cel files
##############################################################################

##############################
#Creat function to find count that are less than 0.05
###############################

count.pscalls <- function(x){
length (which(x<0.05))
}

count_pscalls_threshold <- apply(callps, 1, count.pscalls)
x <- ((which(count_pscalls_threshold>=200))) #-----find probesets that are present in 200 or more with 0.05 P-values
print(length(x)) #516081

probeset_exon_summaries <- rma(rawData, target="probeset", background=TRUE, normalize=TRUE ) #based 

matrixexon_probeset<-exprs(probeset_exon_summaries)
getexon_expression<-as.data.frame(matrixexon_probeset)

print(length(match((row.names(as.data.frame(x))),(row.names(getexon_expression)))))

index_probeset_0.05pval_200<-match((row.names(as.data.frame(x))),(row.names(getexon_expression)))
filt.getexon_expression<-getexon_expression[index_probeset_0.05pval_200,]

##############################################################################
#Print filtered Probesets to file
##############################################################################
write.table(filt.getexon_expression, "PD_RMAnormalized_exon_filtered.txt",col.names=NA)

