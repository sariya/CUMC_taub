#!/bin/Rscript

#
#Sanjeev Sariya
#Date 11/06/2018
#PD Karen Marder, Giuseppe Tosto
#

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
expression_exon<-fread("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/rscripts/PD_RMAnormalized_probe.txt",header=TRUE)

expression_gene<-as.data.frame(expression_gene)
expression_exon<-as.data.frame(expression_exon)

print(dim(expression_gene)) #22011
print(dim(expression_exon)) #1411399

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
print(length(keep_exonprobesets)) #1135900

#Keep exon probeset whose annotation information is found in probsethg19 na30
expression_exon.annot<-expression_exon[match(keep_exonprobesets,expression_exon$V1),]
print(dim(expression_exon.annot)) #1135900

########################################################################
#Filter probeset annotation for filtered genes only
#########################################################################
print(length(intersect(filtered.gene_expression$V1,probesets$transcript_cluster_id))) #17873

#annotation probesets
probesets_transcript_found.geneexpression<-probesets[ which(probesets$transcript_cluster_id %in% gene_to_keep_found_probesetannot.transcript),]

#annotation probesets
print(dim(probesets_transcript_found.geneexpression)) #582790
length(unique(probesets_transcript_found.geneexpression$transcript_cluster_id))#17873

#--exon probesets
print(length(intersect(expression_exon.annot$V1,probesets_transcript_found.geneexpression$probeset_id))) #582790

exonprobesets_genesfound<-intersect(expression_exon.annot$V1,probesets_transcript_found.geneexpression$probeset_id) #582790

exon_probsets_genefounds<-expression_exon.annot[match(exonprobesets_genesfound,expression_exon.annot$V1),]
print(dim(exon_probsets_genefounds)) #582790

#
#For filtered genes and exons now get per exon expression
#

exon_splicing_index<- as.data.frame(matrix(nrow=nrow(exon_probsets_genefounds), ncol=224))

#
#Run loop for each exon filtered 
#Find index of the exon probeset in annotated data 
#Then find transcript of the exon. Find gene expression of the transcript of exon probeset 
#get ratio of exon to the gene expression
#
for(i in 1:nrow(exon_probsets_genefounds)) {

#--get index of probeset from annotation file
index_probeset<- which(probesets_transcript_found.geneexpression$probeset_id== exon_probsets_genefounds[i,1])

if(length(index_probeset) !=0 ){

#--get transcript id of the index found probeset id
transcriptid<-probesets_transcript_found.geneexpression[index_probeset,"transcript_cluster_id"]

#find index of the transcript id found in gene expression data
index_transcript_matched_gene.expression<-which(filtered.gene_expression$V1==transcriptid)

if(length(index_transcript_matched_gene.expression)!=0){

#--calculate ratio of the exon with respect to the gene within sample
exon_splicing_index[i,1:224]<-exon_probsets_genefounds[i,2:225]/filtered.gene_expression[index_transcript_matched_gene.expression,2:225]
}
else{
exon_splicing_index[i,1:224]<-NA
}

}
else{
exon_splicing_index[i,1:224]<-NA
}


}
#--for loop ends
print(dim(exon_splicing_index)) 

row.names(exon_splicing_index)<-exon_probsets_genefounds[,1]
cel_names<-grep("CEL",colnames(joined_annot_exon_filtered),value=TRUE)


