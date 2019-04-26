#!/bin/Rscript

#
#Sanjeev Sariya
#04/24/2019
#Plot using GEMMA GWAS output. The SNPs are harmonized
#

library(Cairo)
library(qqman)

file.name<-"sorted_gwastangles_apoe4"
df.strs<-read.table(file.name,header=TRUE)
print(dim(df.strs))

beta <-df.strs$beta
se <-df.strs$se

lambda = median ((beta/se)^2)/0.4549
lambda

##########################################
##########################################
##########################################

df.strs$rs<-as.character(df.strs$rs)
temp<-strsplit(df.strs$rs,"\\:") #split column for RS
rs_chrmatrix <- do.call(rbind,temp) #do binding
rs_chrmatrix <- as.data.frame(rs_chrmatrix,stringsAsFactors = FALSE)
colnames(rs_chrmatrix)<-c("CHR","BP","SNP")
rs_chrmatrix<-data.frame(rs_chrmatrix,df.strs[,-c(1)])
rs_chrmatrix$CHR<-as.numeric(as.character(rs_chrmatrix$CHR))
rs_chrmatrix$BP<-as.numeric(as.character(rs_chrmatrix$BP))


CairoJPEG(filename = "tangles_qq.jpeg",quality = 75,height=700,width=700)
qq(rs_chrmatrix$p_wald,main="Tangles", xlim = c(0, 7), ylim = c(0,7), pch = 18,col = "blue4", cex = 1.5, las = 1)
dev.off()

