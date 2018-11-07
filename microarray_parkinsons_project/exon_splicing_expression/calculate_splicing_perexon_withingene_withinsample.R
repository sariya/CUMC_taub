#!/bin/Rscript

#
#Sanjeev Sariya
#Date 11/06/2018
#PD Karen Marder, Giuseppe Tosto
#

library(dplyr)
library(purrr)
library(data.table)

transcripts<-fread("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/HuEx-1_0-st-v2.na36.hg19.transcript.csv",header=TRUE,sep=",",skip=23)
transcripts=as.data.frame(transcripts)
print(dim(transcripts))

probesets <- fread("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/exon_level/APT_analysis_scripts_data/HuEx-1_0-st-v2.na30.hg19.probeset.csv",skip=19,header=TRUE,sep=",")
probesets =as.data.frame(probesets)
print(dim(probesets))  #1422046

probesets<-probesets[which(probesets$crosshyb_type==1),]
print(dim(probesets)) #1135900

expression_gene<-fread("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/PD_RMAnormalized_core.txt",header=TRUE)
expression_exon<-fread("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/rscripts/PD_RMAnormalized_exon_filtered.txt",header=TRUE)

expression_gene<-as.data.frame(expression_gene)
expression_exon<-as.data.frame(expression_exon)

print(dim(expression_gene)) #22011
print(dim(expression_exon)) #516081

###########################################################################
#Fix CEL ids from filtered exon expression data. replace - with .
#https://stackoverflow.com/questions/39997273/r-combine-several-gsub-function-ina-pipe
#####################################################################################################

temp_exoncolnames<-colnames(expression_exon)

for(i in 1:length(temp_exoncolnames)){
ifelse((grep("-",temp_exoncolnames[i]) ),temp_exoncolnames[i]<-gsub("-",".",temp_exoncolnames[i]),temp_exoncolnames[i]<-temp_exoncolnames[i] )
}

#
#Filtered Probes don't have X in them. Fix it.
#

for(i in 1:length(temp_exoncolnames)){
ifelse((grep("CEL",temp_exoncolnames[i]) ),temp_exoncolnames[i]<-paste("X",temp_exoncolnames[i],sep=""),temp_exoncolnames[i]<-temp_exoncolnames[i] )
}

colnames(expression_exon)<-temp_exoncolnames
print(length(intersect(expression_gene$V1,transcripts$transcript_cluster_id))) #22011

###############################################################
#Filter gene expression
#---keep genes that are present in probeset cluster ids
###############################################################

print( length(intersect(expression_gene$V1,probesets$transcript_cluster_id))) #17873
gene_to_keep_found_probesetannot.transcript<-(intersect(expression_gene$V1,probesets$transcript_cluster_id))

print(length(match(gene_to_keep_found_probesetannot.transcript,expression_gene$V1))) #17873
filtered.gene_expression<-expression_gene[match(gene_to_keep_found_probesetannot.transcript,expression_gene$V1),]
print(dim(filtered.gene_expression))

#
#Keep exon probes that are found in the probesets crosshyb_type==1
#

#exon probeset means normalized at probeset using RMA
keep_exonprobesets<-intersect(expression_exon$V1,probesets$probeset_id)
print(length(keep_exonprobesets)) #379374

#Keep exon probeset whose annotation information is found in probsethg19 na30
expression_exon.annot<-expression_exon[match(keep_exonprobesets,expression_exon$V1),]
print(dim(expression_exon.annot)) #379374

########################################################################
#Filter probeset annotation for filtered genes only
#########################################################################
print(length(intersect(filtered.gene_expression$V1,probesets$transcript_cluster_id))) #17873

#annotation probesets
probesets_transcript_found.geneexpression<-probesets[ which(probesets$transcript_cluster_id %in% gene_to_keep_found_probesetannot.transcript),]

print(dim(probesets_transcript_found.geneexpression)) ##582790     39

#annotation probesets

print(length(unique(probesets_transcript_found.geneexpression$transcript_cluster_id)))#17873

##############################
#--exon probesets
##########################
#--find the ones for which annotation is present

print(length(intersect(expression_exon.annot$V1,probesets_transcript_found.geneexpression$probeset_id))) #254798
exonprobesets_genesfound<-intersect(expression_exon.annot$V1,probesets_transcript_found.geneexpression$probeset_id) #

exon_probsets_genefounds<-expression_exon.annot[match(exonprobesets_genesfound,expression_exon.annot$V1),]
print(dim(exon_probsets_genefounds)) #254798

##
#Filter probeset annotation based on exon probsets availalbe
##

print(length(intersect(exon_probsets_genefounds$V1,probesets_transcript_found.geneexpression$probeset_id))) #254798
probsets_filtered.annotation_foundexonexpression<-intersect(exon_probsets_genefounds$V1,probesets_transcript_found.geneexpression$probeset_id)
print(length(probsets_filtered.annotation_foundexonexpression))

filtered_annotation_foundexonexpression<-probesets_transcript_found.geneexpression[match(probsets_filtered.annotation_foundexonexpression,probesets_transcript_found.geneexpression$probeset_id),]
print(dim(filtered_annotation_foundexonexpression))

genes_further_filtered<-intersect(filtered.gene_expression[,1],filtered_annotation_foundexonexpression$transcript_cluster_id)
filtered.gene_expression2<-filtered.gene_expression[match(genes_further_filtered,filtered.gene_expression[,1]),]
#
#Update colnmae of gene expressions to merge and later get exon ratio
#

newcolnames_geneexpression<-c("geneExpression_transcriptID", colnames(filtered.gene_expression2)[grep("CEL",colnames(filtered.gene_expression2))])
colnames(filtered.gene_expression2)<-newcolnames_geneexpression # 16941   225

#
#For filtered genes and exons now get per exon expression
#

###https://stackoverflow.com/questions/28765268/how-to-do-iterative-joins-in-plyr##

join_exon.geneExprs.probeids_transcriptids<-(left_join(exon_probsets_genefounds,filtered_annotation_foundexonexpression,by=c("V1"="probeset_id")) %>% left_join(.,filtered.gene_expression2,by=c("transcript_cluster_id"="geneExpression_transcriptID")))

print(length(grep("CEL",colnames(join_exon.geneExprs.probeids_transcriptids))))

##https://stackoverflow.com/questions/41452097/match-by-id-and-divide-column-values-across-two-dataframes ##

exonexp.ratio_with_geneexpres<-cbind(id=join_exon.geneExprs.probeids_transcriptids[,1], map2_df(join_exon.geneExprs.probeids_transcriptids[,2:225], join_exon.geneExprs.probeids_transcriptids[,264:487], `/`))

#
#For filtered genes and exons now get per exon expression
#

print(dim(exonexp.ratio_with_geneexpres))  ##254798
print("we have calculated ratio for SI")

#
#Print spliced exon ratio to files
#make one chunk of 25000 exon
#

removed.x_colnames<-colnames(exonexp.ratio_with_geneexpres) %>% gsub(".x","",.)
colnames(exonexp.ratio_with_geneexpres)<-removed.x_colnames

list_split_indexExon<-split(exonexp.ratio_with_geneexpres, (seq(nrow(exonexp.ratio_with_geneexpres))-1) %/% 25000) 
length(list_split_indexExon) # 11 #

for(i in 1:length(list_split_indexExon) )
{
file<-paste("exon_chunk_ratio",i,sep="_") #--make file name for chunk
file_dest<-paste("./ratio_exon_to_gene",file,sep="/")
write.table(list_split_indexExon[[i]], file_dest, sep="\t", quote=F, row.names=FALSE,col.names=TRUE)
}




"R version 3.4.2 (2017-09-28)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 9 (stretch)

Matrix products: default
BLAS: /mnt/mfs/cluster/bin/R-3.4/lib/libRblas.so
LAPACK: /mnt/mfs/cluster/bin/R-3.4/lib/libRlapack.so

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] data.table_1.11.8 purrr_0.2.5       dplyr_0.7.7

loaded via a namespace (and not attached):
 [1] tidyselect_0.2.5 compiler_3.4.2   magrittr_1.5     assertthat_0.2.0
 [5] R6_2.3.0         pillar_1.3.0     bindrcpp_0.2.2   glue_1.3.0
 [9] tibble_1.4.2     crayon_1.3.4     Rcpp_0.12.19     pkgconfig_2.0.2
[13] rlang_0.3.0.1    bindr_0.1.1
"


#SD_perprobeset<-as.data.frame(apply(exon_probsets_genefounds[,2:225],1,sd))
#SD_perprobeset$probeset_id<-exon_probsets_genefounds$V1
#colnames(SD_perprobeset)<-c("sd","V1")

#sd_perprobeset_joined_annot.transcript<-left_join(SD_perprobeset,filtered_annotation_foundexonexpression,by=c("V1"="probeset_id"))
#median_pertranscript_fromprobeset.sd<-as.data.frame(sd_perprobeset_joined_annot.transcript %>% group_by(transcript_cluster_id) %>% summarize_at(vars(sd),funs(median))  %>% 
#rename(median_value=sd))  %>% 
#left_join(.,sd_perprobeset_joined_annot.transcript,by=c("transcript_cluster_id"="transcript_cluster_id"))
#
#Normalized SI
#
#joined_SI_median_withintranscript<-left_join(exonexp.ratio_with_geneexpres,median_pertranscript_fromprobeset.sd,by=c("id"="V1"))
###########