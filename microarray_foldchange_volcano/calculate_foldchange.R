#!/bin/Rscript

#
#Date 11/06/2018
#sanjeev Sariya
#Use amanda myers data to learn fold change,. design and other volcance plots
#Do fold change and volcane plot
#
#
library(dplyr)

#--read matrix already adjusted and regressed
summed<-read.table("duplicate_genes_summed",header=TRUE)
print(dim(summed))

#get gene names
gene_names<-summed[,1]

row.names(summed)<-gene_names

#--make group, N are normal, A ends for AD 
group<-factor(c(rep(0,length(grep("N$",colnames(summed)))),rep(1,length(grep("A$",colnames(summed)))) ))
group <- relevel(group, "0")

design<-model.matrix(~group+0)

fit <- lmFit(summed, design)
fit2 <- eBayes(fit)
topTable(fit2, coef = 2)

results <- decideTests(fit2)
a <- vennCounts(results)
#mfrow.old <- par()$mfrow
#par(mfrow=c(1,2))
#vennDiagram(a)
#vennDiagram(results, include=c("up", "down"), counts.col=c("red", "blue"), circle.col = c("red", "blue", "green3")) #par(mfrow=mfrow.old)
