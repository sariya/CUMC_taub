#!/bin/Rscript

#
#Date 12/17/2018
#Sanjeev Sariya
#Dr. Giuseppe Tosto, Karen Marder
#
######################################
# Gene level normalized data
######################################

## https://cran.r-project.org/web/packages/myTAI/vignettes/Expression.html
#
#Read pheno data
#Get case and controls calculate mean per gene in case, then calculate mean per gene in control

#R 3.4.2
library(data.table) #data.table 1.11.8 
library(dplyr) #dplyr_0.7.7
library(huex10sttranscriptcluster.db) # huex10sttranscriptcluster.db_8.7.0 org.Hs.eg.db_3.5.0


file.expression<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/gene_expr_rin_site_adjusted/PD_RMAnormalized_core.txt"
file.pheno<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/gene_expr_rin_site_adjusted/pheno_rin_study"

pd<-data.table::fread(file.expression, showProgress = TRUE) 
transposed<-t(pd)

gene_ids<-pd$V1 #get gene ids
colnames(transposed)<-gene_ids

pheno<-read.table( file.pheno ,header=TRUE)
row.names(pheno)<-pheno$IID
PDmerge<-merge(transposed,pheno, by=0,all=TRUE)
#find missing Cel pheno IIDs
missing_celpehno<-which(is.na(PDmerge$IID))
PDmerge<-PDmerge[-c(missing_celpehno),] #remove data for missing cel phenotypes

PDmerge<-PDmerge[,-c(1)]

print(dim(PDmerge))
PDmerge$AGE_SCALE<-(scale(PDmerge$AGE))[,1]
names<-colnames(PDmerge)

#
# https://stackoverflow.com/questions/9723208/aggregate-summarize-multiple-variables-per-group-e-g-sum-mean
#

#calculate mean of each gene for case and contol
mean_case_control<-as.data.frame(PDmerge %>% group_by(PD) %>% summarise_at(vars(-AGE_SCALE, -IID, -CATEGORY, -SEX, -AGE, -PD, -LRKK2, -RIN, -FID, -STUDY), mean ))
print(dim(mean_case_control))

mean_case_control[3,]<-NA
mean_case_control[3,2:ncol(mean_case_control)]<-mean_case_control[1,2:ncol(mean_case_control) ]- mean_case_control[2,2:ncol(mean_case_control) ]
rownames(mean_case_control)<-c("Control","Case","logFC")

transposed_logfc<-t(mean_case_control)
transposed_logfc<-transposed_logfc[-c(1),]

transposed_logfc<-as.data.frame(transposed_logfc)

#https://www.statmethods.net/management/sorting.html
sorted_transposed_logfc<- transposed_logfc[order(transposed_logfc$logFC),] #sort df
print(dim(sorted_transposed_logfc))

sorted_transposed_logfc$PROBEID<-rownames(sorted_transposed_logfc) #use this in joining later

############################################################################
#Get gene names 
############################################################################

collapser <- function(x){
x %>% unique %>% sort %>% paste(collapse = "|")
}  

#https://bioconductor.org/packages/devel/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf

probe.annots  <- as.data.frame(AnnotationDbi::select(
x       = huex10sttranscriptcluster.db,
keys    = as.character( rownames(sorted_transposed_logfc)),
columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
keytype = "PROBEID"  ) %>%
group_by(PROBEID) %>% summarise_each(funs(collapser)) %>%   ungroup)
  
print(dim(probe.annots))
probe.annots  <-probe.annots[which(probe.annots$SYMBOL!=""  ),] 

# drop columns https://stackoverflow.com/questions/5234117/how-to-drop-columns-by-name-in-a-data-frame
genenames_foldchange<-subset( dplyr::left_join(sorted_transposed_logfc ,probe.annots,by = c("PROBEID")) , select=-c(ENTREZID)) 
genenames_foldchange_noNAs<-genenames_foldchange[complete.cases(genenames_foldchange),]
write.table(genenames_foldchange_noNAs,file="geneexpression.sorted_foldchange.tsv",sep="\t",append = FALSE, quote = FALSE, row.names = FALSE,col.names = TRUE)


#get mean age per group
as.data.frame(PDmerge %>% group_by(PD) %>% summarise_at(vars(AGE), funs(mean(., na.rm = TRUE)) ) ) 
#  PD      AGE
#  0 58.06731
#  1 69.04274

#sd of age in case and control
as.data.frame(PDmerge %>% group_by(PD) %>% summarise_at(vars(AGE), funs(sd(., na.rm = TRUE)) ) )
#  PD       AGE
#  0 17.862617
#  1  9.695709


