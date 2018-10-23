#!/bin/Rscript

#
#Date 10/22/2018
#
#Sanjeev Sariya
#
#heatmap for PD data. Pick Pvalue based on your need
#

library(dplyr)
library(Cairo)
library(gplots)
library(RColorBrewer)
library(marray)
library(gtools)

#Source of code mainly
##http://www.hiv.lanl.gov/content/sequence/HEATMAP/heatmap.html?sample_input=1
#
pheno<-read.table( "/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/exon_level/APT_analysis_scripts_data/OUT_split_exon_modeling/pheno_rin_study" ,header=TRUE)
row.names(pheno)<-pheno$IID

file.expression<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/PD_RMAnormalized_core.txt"
file.model_genes<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/gene_expr_rin_site_adjusted/model3_genes"

df.expression<-data.table::fread(file.expression, showProgress = TRUE,header=TRUE)
print(dim(df.expression))

df.model.genes<-data.table::fread(file.model_genes, showProgress = TRUE,header=TRUE)
print(dim(df.model.genes))

pvalue_threshold<-0.005 #suit based on your need - model 1 and model 2 are good with this
print(length(which(df.model.genes$Pvalue<pvalue_threshold)))

genes_pvalue_threshold<- df.model.genes[which(df.model.genes$Pvalue<pvalue_threshold),]
print(dim(genes_pvalue_threshold))

#--split wierd /\/\ things to extract gene names
split_gene_assignment<-do.call (rbind,strsplit( unlist(genes_pvalue_threshold[,10]) ,"\\/\\/"))
split_gene_assignment<-split_gene_assignment[,2]

genes_pvalue_threshold$gene_names<-unlist(lapply(split_gene_assignment, trimws))

print(length(genes_pvalue_threshold$gene_names))
print(length(unique(genes_pvalue_threshold$gene_names)))
print(dim(genes_pvalue_threshold))

genes_pvalue_threshold_genes<-genes_pvalue_threshold[,c(1,2,20)]
which(genes_pvalue_threshold_genes$gene_names=="---")
print(dim(genes_pvalue_threshold_genes))

#get rid of genes with names as ---
genes_pvalue_threshold_genes<-genes_pvalue_threshold_genes[-c(which(genes_pvalue_threshold_genes$gene_names=="---")),]

#--check for duplicate gene names
print(length(genes_pvalue_threshold_genes$gene_names))
print(length(unique(genes_pvalue_threshold_genes$gene_names)))

print(dim(genes_pvalue_threshold_genes))

#do a join on porbe ids
expression_value.shortlisted_genes<-left_join(genes_pvalue_threshold_genes,df.expression,by=c("GeneProbeset"="V1"))
row.names(expression_value.shortlisted_genes)<-expression_value.shortlisted_genes$gene_names

cels_genes<-expression_value.shortlisted_genes[,-c(1,2,3)]
print(dim(cels_genes))
tranposed_celsgenes<-t(as.matrix(cels_genes))
print(dim(tranposed_celsgenes))

cels_merge<-merge(tranposed_celsgenes,pheno, by=0,all=TRUE)
print(dim(cels_merge))

#--get rid of poor IIDs
missing_celpheno<-which(is.na(cels_merge$IID))
cels_merge<-cels_merge[-c(missing_celpheno),]

##control - green - cases red
color.map <- function(PD) { 
if (PD==0) "#ADFF2F"  else "#FF0000"  
}

#--use patient colors for side coloring of columns
patientcolors <- lapply(cels_merge$PD, color.map)

rownames(cels_merge)<-cels_merge[,1]
print(dim(cels_merge))

#--get rid of columns not needed from pheno merging
cels_merge<-cels_merge[ , -which(names(cels_merge) %in% c("Row.names","IID","CATEGORY","SEX","AGE","PD","LRKK2","RIN","FID","STUDY"))]
cels_merge_trans<-as.matrix(t(cels_merge))
print(dim(cels_merge))

#--set color key
vec_m <- c(cels_merge_trans)

#--make breaks. and use them to get color legend
brks = seq(min(vec_m, na.rm=T), max(vec_m, na.rm=T), by=diff(range(vec_m, na.rm=T))/9)

#--set color counts
colors_count = length(brks) - 1

#--set your color key higher the darker

my_palette =brewer.pal(colors_count,"YlOrRd")

CairoPNG("model3_16genes.png",height=1200,width=1200)

heatmap.2((cels_merge_trans),
key.xlab="Log gene expression", key.ylab="Count",
ColSideColors=unlist(patientcolors),
labCol = FALSE,breaks=brks,   key=TRUE,trace="none", 
distfun=function(x) dist(x,method="euclidean") , hclustfun=function(m) hclust(m,method="ward.D"), margins=c(5,12),
 xlab = "Samples",ylab = "Genes"	 ,cexRow=1.5 ,
col=my_palette ,
 main="Top Genes Model 3"
)

legend("topright",      
    legend = c("Control","Case"),
    col = c( "#ADFF2F","#FF0000"  ), 
    lty= 1, lwd = 5,    cex=1.2
    )

dev.off()
######################

#
"
R version 3.4.2 (2017-09-28)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 9 (stretch)

Matrix products: default
BLAS: /mnt/mfs/cluster/bin/R-3.4/lib/libRblas.so
LAPACK: /mnt/mfs/cluster/bin/R-3.4/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] gtools_3.8.1       marray_1.58.0      limma_3.34.9       RColorBrewer_1.1-2
[5] gplots_3.0.1       Cairo_1.5-9        dplyr_0.7.7

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.19       crayon_1.3.4       assertthat_0.2.0   bitops_1.0-6
 [5] R6_2.3.0           magrittr_1.5       KernSmooth_2.23-15 pillar_1.3.0
 [9] rlang_0.3.0        gdata_2.18.0       bindrcpp_0.2.2     glue_1.3.0
[13] purrr_0.2.5        compiler_3.4.2     pkgconfig_2.0.2    caTools_1.17.1.1
[17] bindr_0.1.1        tidyselect_0.2.5   tibble_1.4.2

"


