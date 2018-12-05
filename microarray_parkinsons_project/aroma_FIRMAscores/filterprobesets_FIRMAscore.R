#!/bin/Rscript

#
#Date 11/15/2018
#
#Sanjeev Sariya
#

probesets <- fread("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/ratio_exon_to_gene/HuEx-1_0-st-v2.na30.hg19.probeset.csv",skip=19,header=TRUE,sep=",")
probesets <-as.data.frame(probesets)
print(dim(probesets))  #1422046

probesets<-probesets[(which(probesets$crosshyb_type==1)),]
dim(aromaprobesets)

print(length(intersect(aromaprobesets$groupName,probesets$probeset_id)))

matched_probesets<-intersect(aromaprobesets$groupName,probesets$probeset_id)
index_probesets<-match(matched_probesets,aromaprobesets$groupName)
aromaprobesets.crosshyb<-aromaprobesets[index_probesets,]

print(dim(aromaprobesets.crosshyb))

write.table(aromaprobesets.crosshyb,"firmascores_colnamesfixed.filteredcross",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

list_split_firmaoutput<-split(aromaprobesets.crosshyb, (seq(nrow(aromaprobesets.crosshyb))-1) %/% 30000) 
print(length(list_split_firmaoutput))

#
#Split into multiple chunks
#

for(i in 1:length(list_split_firmaoutput) )
{
file<-paste("aroma.firmascore_chunk",i,sep="_") #--make file name for chunk
file_dest<-paste("./",file,sep="/")
write.table(list_split_firmaoutput[[i]], file_dest, sep="\t", quote=F, row.names=FALSE,col.names=TRUE)
}



