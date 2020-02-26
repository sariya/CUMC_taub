#!/usr/bin/env Rscript	

#
#Date 02/26/2019
#Sanjeev Sariya
#
## /mnt/mfs/cluster/bin/R-3.4/bin/Rscript /home/ss5505/scripts_cumc/automate_analyses/merge_gmmatresults_withgrch37_genes/merge_gmmat_withgenes.R -r merged.sorted -o model1_hgwas123467PR_genes.sorted
##
suppressMessages(library(dplyr))
suppressMessages(library(argparse))

parser <- ArgumentParser(description="merge gmmat output with gene-chr-start-position")

parser$add_argument('-r',"--gmmat",help="File with pvalue from GMMAT analyses",required=TRUE)
parser$add_argument('-o',"--output",help="name of output file",required=TRUE)
args <- parser$parse_args() #make it a data structure

file.results<-normalizePath(args$gmmat)
file.output<-args$output

file.genes<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/scripts/ANNOTATION/PLINK_annotation/NCBI_37_v3_GRCh37_forPLINK"

df.genes<-read.table(file.genes,header=FALSE)
df.results<-read.table(file.results,header=TRUE)

print(dim(df.genes))
print(dim(df.results))
colnames(df.genes)<-c("Chr","Start","End","GeneName")

lefted.join<-left_join(df.results,df.genes,by=c("group"="GeneName"))

lefted.join <-subset(lefted.join,select=c("group", "Chr", "Start","End","n.variants", "miss.min","miss.mean",
 "miss.max","freq.min", "freq.mean", "freq.max","B.score","B.var",
"B.pval","S.pval","O.pval", "O.minp", "O.minp.rho",
"E.pval"))

colnames(lefted.join )<-c("Gene", "Chr", "Start","End","n.variants", "miss.min","miss.mean",
"miss.max","freq.min", "freq.mean", "freq.max","B.score","B.var",
"B.pval","S.pval","O.pval", "O.minp", "O.minp.rho",
"E.pval")

lefted.join <-lefted.join [order(lefted.join$O.pval ),]

write.table(x=lefted.join , file = file.output, append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE)





