#!/usr/bin/env Rscript

###
#Date 03/11/2020
###Sanjeev Sariya
#####################################
##
## "SNP","CHR", "START","END","P" these columns are needed in file
## R 3.4
## Rscript create_manhattan_genebased.R -f plot_genes.input -x 6 -t "Gene based: Rare SNPs Model 1HGWAS123467 PR" -o manh_plot_rareSNPs_model1_HGWAS1234567PR
###  awk '{if(FNR==1){print "SNP","CHR", "START","END","P"} else{print  $1,$2, $3,$4,$16 }}' sortedgenenames_positions_rare_model1_merged.txt > plot_genes.input
#####################################

suppressMessages(library(qqman))
library(Cairo)
library(argparse)

parser <- ArgumentParser(description="Create Manh plot using parameters")

parser$add_argument('-f',"--file",help="File with pvalue and other information",required=TRUE)
parser$add_argument('-x',"--axis",help="max axis label 5, 6, 7",required=TRUE)
parser$add_argument('-o',"--outname",help="output name of plot",required=TRUE)
parser$add_argument('-t',"--title",help="header you'd like to see in plot",required=TRUE)

args <- parser$parse_args() #make it a data structure

file_genebased<-normalizePath(args$file) ##file with pvalue
axis_limit<-as.numeric(args$axis) ##column number that contains pvalue
title<- args$title ##the title you'd like 
output_filename<-as.character(args$outname)

df.genepvalue<-read.table(file_genebased,header=TRUE)
print(dim(df.genepvalue)) ##check if data are loaded fine ## SNP CHR START END PVAL

df.genepvalue$BP<-(df.genepvalue$START+df.genepvalue$END)/2

#Use qqman to plot manht graph
CairoPNG(filename = paste(output_filename,".png",sep=""),quality = 75, height=800, width=900)

par(mar = c(5,6,4,1))

suppressWarnings(manhattan(df.genepvalue,main= title, 
ylim = c(0,axis_limit),  cex = 1.5, cex.main=1.5,
cex.axis = 1.2, cex.lab = 1.7,font.lab=2,
col = c("goldenrod", "gray0"),
genomewideline = -log10(2.2E-06),
suggestiveline = F,
chrlabs = c(1:22)))
dev.off()


##col = c("#A759B1", "#180B5D", "#660E72", "#7468BA", "#5B2392"),