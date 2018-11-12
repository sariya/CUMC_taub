#!/bin/Rscript

#
#Date 11/12/2018
#
#Sanjeev Sariya
#
#heatmap for PD data. Pick Pvalue based on your need
#exon level data

library(huex10stprobeset.db)
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

file.expression<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/rscripts/allexons_gene_ratio.txt"
file.model_exon<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/ratio_exon_to_gene/RIN_adjusted/Model4/adjusted_valuesOverall"

df.expression<-data.table::fread(file.expression, showProgress = TRUE,header=TRUE)
print(dim(df.expression))

df.model.exons<-data.table::fread(file.model_exon, showProgress = TRUE,header=TRUE)
print(dim(df.model.exons))

#temp_exoncolnames<-colnames(df.expression)

###########################################################################
#Fix CEL ids from filtered exon expression data. replace - with .
#https://stackoverflow.com/questions/39997273/r-combine-several-gsub-function-ina-pipe
#####################################################################################################

##for(i in 1:length(temp_exoncolnames)){ ifelse((grep("-",temp_exoncolnames[i]) ),temp_exoncolnames[i]<-gsub("-",".",temp_exoncolnames[i]),temp_exoncolnames[i]<-temp_exoncolnames[i] ) }
#
#Filtered Probes don't have X in them. Fix it.
#
##for(i in 1:length(temp_exoncolnames)){ ifelse((grep("CEL",temp_exoncolnames[i]) ),temp_exoncolnames[i]<-paste("X",temp_exoncolnames[i],sep=""),temp_exoncolnames[i]<-temp_exoncolnames[i] ) }
#colnames(df.expression)<-temp_exoncolnames

print(dim(df.model.exons))
##[1] 254798     10

#get adjusted ones
df.model.exons<-df.model.exons[(which(df.model.exons$adjusted_Overall_pvalue <0.05)),]
print(dim(df.model.exons))
##[1] 1908   10

print(min(df.model.exons$adjusted_cat4)) #0.002737049
print(min(df.model.exons$adjusted_cat3)) #4.797503e-05
print(min(df.model.exons$adjusted_cat2)) #4.909713e-11

df.model.exons<-df.model.exons[(which(df.model.exons$adjusted_cat2<0.05)),]
print(dim(df.model.exons)) ## 1232   10

collapser <- function(x){
x %>% unique %>% sort %>% paste(collapse = "|")
}  

probe.annots  <- as.data.frame(AnnotationDbi::select(
x       = huex10stprobeset.db,
keys    = as.character(df.model.exons$Probeset),
columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
keytype = "PROBEID"  ) %>%
group_by(PROBEID) %>% summarise_each(funs(collapser)) %>%   ungroup)

#--delete blank rows  
probe.annots  <-probe.annots[which(probe.annots$SYMBOL!="" & probe.annots$ENSEMBL!="") ,]
 
#--delete genes
probe.annots.z <- probe.annots[order(probe.annots$SYMBOL),] 
probe.annots.z <- probe.annots.z[!duplicated(probe.annots.z$SYMBOL),]

probe.annots.z$PROBEID<-as.numeric(probe.annots.z$PROBEID)
joinedpvalue_adjusted_probe.annot<-left_join(df.model.exons,probe.annots.z,by=c("Probeset"="PROBEID"))
joinedpvalue_adjusted_probe.annot<-joinedpvalue_adjusted_probe.annot[complete.cases(joinedpvalue_adjusted_probe.annot),]

write.table(joinedpvalue_adjusted_probe.annot,"pvalues.adjusted.geneSymbols_EnsemblIDs",quote=FALSE,row.names=FALSE,col.names=TRUE)

print(dim(joinedpvalue_adjusted_probe.annot))
## 516  13

#--keep only limited

joinedpvalue_adjusted_probe.annot<-joinedpvalue_adjusted_probe.annot[ , which(colnames(joinedpvalue_adjusted_probe.annot) %in% 
c("Probeset","ENSEMBL","SYMBOL"))]

#do a join on probe ids
expression_value.shortlisted_genes<-left_join(joinedpvalue_adjusted_probe.annot,df.expression,by=c("Probeset"="id"))
row.names(expression_value.shortlisted_genes)<-expression_value.shortlisted_genes$SYMBOL

#--get rid of  Probeset         ENSEMBL  SYMBOL
cels_genes<-expression_value.shortlisted_genes[ , -which(colnames(expression_value.shortlisted_genes) %in% c("Probeset","SYMBOL","ENSEMBL"))]
print(dim(cels_genes))
tranposed_celsgenes<-t(as.matrix(cels_genes))
print(dim(tranposed_celsgenes))

cels_merge<-merge(tranposed_celsgenes,pheno, by=0,all=TRUE)
print(dim(cels_merge))

#--get rid of poor IIDs
missing_celpheno<-which(is.na(cels_merge$IID))
cels_merge<-cels_merge[-c(missing_celpheno),]
print(dim(cels_merge))

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

CairoPNG("Poissonmodel4_topcat2_pval_adjusted2.png",height=1200,width=1200)
heatmap.2((cels_merge_trans),cexRow=1.2,
key.xlab="Ratio of exon to gene expression", key.ylab="Count",ColSideColors=unlist(patientcolors),labCol = FALSE,breaks=brks,key=TRUE,trace="none", 
distfun=function(x) dist(x,method="euclidean") , hclustfun=function(m) hclust(m,method="complete"), margins=c(5,21),
xlab = "Samples",ylab = "Genes",col=my_palette ,
main="Top Genes Model 4 Category 2")

legend("topright", legend = c("Control","Case"),
col = c( "#ADFF2F","#FF0000"  ), lty= 1, lwd = 5,cex=1.2)

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
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] gtools_3.8.1                       marray_1.58.0
 [3] limma_3.34.9                       RColorBrewer_1.1-2
 [5] gplots_3.0.1                       Cairo_1.5-9
 [7] dplyr_0.7.7                        huex10sttranscriptcluster.db_8.7.0
 [9] org.Hs.eg.db_3.5.0                 AnnotationDbi_1.42.1
[11] IRanges_2.12.0                     S4Vectors_0.16.0
[13] Biobase_2.38.0                     BiocGenerics_0.24.0

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.19       compiler_3.4.2     pillar_1.3.0       bindr_0.1.1
 [5] bitops_1.0-6       digest_0.6.18      bit_1.1-14         RSQLite_2.1.1
 [9] memoise_1.1.0      tibble_1.4.2       pkgconfig_2.0.2    rlang_0.3.0.1
[13] DBI_1.0.0          bindrcpp_0.2.2     caTools_1.17.1.1   bit64_0.9-7
[17] tidyselect_0.2.5   data.table_1.11.8  glue_1.3.0         R6_2.3.0
[21] gdata_2.18.0       purrr_0.2.5        blob_1.1.1         magrittr_1.5
[25] assertthat_0.2.0   KernSmooth_2.23-15 crayon_1.3.4

"


