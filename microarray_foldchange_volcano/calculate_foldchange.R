#!/bin/Rscript

#
#Date 11/06/2018
#sanjeev Sariya
#Use amanda myers data to learn fold change,. design and other volcance plots
#Do fold change and volcane plot
#
#
library(dplyr)
library(limma)
library(Cairo)
#--read matrix already adjusted and regressed
summed<-read.table("duplicate_genes_summed",header=TRUE)
print(dim(summed))

#get gene names
gene_names<-summed[,1]

summed<-summed[,-c(1)]
row.names(summed)<-gene_names

#--make group, N are normal, A ends for AD 
group<-factor(c(rep(0,length(grep("N$",colnames(summed)))),rep(1,length(grep("A$",colnames(summed)))) ))
#group <- relevel(group, "0")

design<-model.matrix(~group+0)

fit <- lmFit(summed, design)
fit2 <- eBayes(fit)
topTable(fit2, coef = 2)

results <- decideTests(fit2)
a <- vennCounts(results)

CairoPDF("volcanoplot.pdf")
volcanoplot(fit2 , coef=2, xlim =c(-1,1))
dev.off()