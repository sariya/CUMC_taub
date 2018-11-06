#!/bin/Rscript


setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/CEL/")
####Normalization based on Core 
###Parkinon's expression data
###Date: March  20  2018
###Sanjeev Sariya
###Karen-Giu project

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
geneSummaries <- rma(rawData, target="core", background=TRUE, normalize=TRUE ) 

PDmatrixgenes<-exprs(geneSummaries)
getexpression<-as.data.frame(as.ffdf(PDmatrixgenes))
print("Get exppression complete")

#--leave first column blank
write.table(getexpression, "PD_RMAnormalized_core.txt",col.names=NA)

