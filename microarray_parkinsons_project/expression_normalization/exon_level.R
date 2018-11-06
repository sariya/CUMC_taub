#!/bin/Rscript

setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/CEL/")
####Normalization based on probe  
###Parkinon's expression data
###Date: April 04 2018
###Sanjeev Sariya
###Karen-Giu project
#
#https://www.biostars.org/p/271379/#307326
#
###--------------------------------
###Load libraries
###----------------------------------
library("Cairo")
library("lme4",lib.loc="~/R_LIB/")
library("oligo",lib.loc="~/R_LIB")
library("ff",lib.loc="~/R_LIB")
library("pd.huex.1.0.st.v2",lib.loc="~/R_LIB") 

setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/CEL")
getwd()
celFiles <- list.celfiles("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/CEL")

wafer <- substr(celFiles, 1, 3)
experiment <- substr(celFiles, 5, 5)
tmp <- substr(celFiles, 4, 4)
complex <- rep('+', length(tmp))

info <- data.frame(wafer=wafer, experiment=experiment, complex=complex)
rownames(info) <- celFiles
metadata <- data.frame(labelDescription=c('wafer', 'experiment', 'complex'), channel=factor('_ALL_'))
rawData <- read.celfiles(celFiles)

sampleNames(rawData) <- celFiles
pd <- new('AnnotatedDataFrame', data=info, varMetadata=metadata)
phenoData(rawData) <- pd

###########################################################################
##############normalization method#########################################
###########################################################################
probeset_exon_summaries <- rma(rawData, target="probeset", background=TRUE, normalize=TRUE ) #based 

matrixexon_probeset<-exprs(probeset_exon_summaries)
getexpression_probsetexon<-as.data.frame(as.ffdf(matrixexon_probeset))
print("Get exppression complete")

#--leave first column blank
write.table(getexpression_probsetexon, "PD_RMAnormalized_exon.txt",col.names=NA)

