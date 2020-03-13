#!/usr/bin/env Rscript

#
#Date 02/26/2020
#Sanjeev Sariya
#
#
#Merge data from GEMMA output with RS ids, CHR, position
#
# /mnt/mfs/cluster/bin/R-3.4/bin/Rscript create_qqplot.R -f filename -c column_with_pvalue \
# -x axis_max_value -n "name_of_plot"
# -t "tite_in_plot"
##

### /mnt/mfs/cluster/bin/R-3.4/bin/Rscript  /mnt/mfs/hgrcgrid/homes/ss5505/scripts_cumc/automate_analyses/create_qqplot/create_qqplot_withlambda.R -f merged.sorted -c 13 -o "qq_HGWAS123467PR_model1_commonrareSNPs" -x 5 -t "HGWAS123467PR Model1 Common Rare SNPs"
##/mnt/mfs/cluster/bin/R-3.4/bin/Rscript  /mnt/mfs/hgrcgrid/homes/ss5505/scripts_cumc/automate_analyses/create_qqplot/create_qqplot_withlambda.R  -f gwas_model1_CHdbgenotyped_sexagebraak_disease.tab -c 5 -o "qq_GWASmodel1_CH_sex_age_PCs_braak_disease" -x 4 -t  "GWASmodel1 CH Database with sex age PCs braak disease"
#

library(qqman)
library(ggplot2)
library(Cairo)
library(argparse)

############################################################################
############################################################################
############################################################################

parser <- ArgumentParser(description="Create Q-Q plot using parameters")

parser$add_argument('-f',"--file",help="File with pvalue and other information",required=TRUE)
parser$add_argument('-c',"--col",help="Scol number for P-value",required=TRUE)
parser$add_argument('-x',"--axis",help="max axis label 5, 6, 7",required=TRUE)
parser$add_argument('-o',"--outname",help="name of plot",required=TRUE) ##output file name of PNG file
parser$add_argument('-t',"--title",help="header you'd like to see in plot",required=TRUE)

args <- parser$parse_args() #make it a data structure

file_pvalue<-normalizePath(args$file) ##file with pvalue
column_number<-as.numeric(args$col) ##column number that contains pvalue
axis_limit<-as.numeric(args$axis) ##column number that contains pvalue
title<-args$title ##the title you'd like 
output_filename<-args$outname

############################################################################
############################################################################
############################################################################

##add red color to the conf interval
gg_qqplot <- function(ps, axis_limit.plot, ci = 0.95) {

  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))  )

  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) + geom_ribbon(fill = "grey3",
      mapping = aes(x = expected, ymin = clower, ymax = cupper), alpha = 0.1 ) +
	  ###EA8A4F
    geom_point(aes(expected, observed), shape = 1, size = 1,colour=("navyblue")) + 
    xlim(0,axis_limit.plot) +ylim(0,axis_limit.plot)+
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, size=1.5 ) +
    xlab(log10Pe) + ylab(log10Po)
}


##
#Function ends
##
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

############################################################################
############################################################################
############################################################################

df.pvalue<-read.table(file_pvalue,header=TRUE)
print(dim(df.pvalue)) ##check if data are loaded fine

lambda_aesthetics<-inflation(df.pvalue[,column_number])

image_name <- paste(output_filename,".png",sep="")
CairoPNG(filename = image_name ,quality = 75, height=800, width=900) 

gg_qqplot(df.pvalue[,column_number],axis_limit ) +
theme_bw(base_size = 14) + 
theme( axis.ticks = element_line(size = 0.9),
panel.grid = element_blank(),axis.text.x = 
element_text(face="bold"),axis.text.y = 
element_text(face="bold") ) +  labs(title= title) + 
annotate( geom = "text", x = -Inf,  y = Inf, hjust = -0.15,vjust = 1 + 0.15 * 3, 
label = as.expression(bquote( lambda~"="~.(sprintf("%.2f",lambda_aesthetics) ) ) ), size = 6 )

dev.off()

print("Check current working directory")

