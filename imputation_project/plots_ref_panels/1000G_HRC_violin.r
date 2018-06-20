#!/bin/Rscript

##
#Date 05/16/2018
#
#MAF box plot for KGP and HRC
#
#Read KGP data first. 
#
library(dplyr)
library("data.table")
library("Cairo")
library(plyr)
library(ggplot2)

setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/HRC_WellcomeTrust/imputation_1000indiv/impute2/1KGP/unharmonized/merged_info_files/")
filename<-"merged_info_CHR"
ext<-".info"
df_all <- as.list(rep("", 22)) 

for(i in seq(1,22)){
file_impute2_info<-paste(filename,i,ext,sep="")

newdf<-paste("df",i,sep="")

df_all[[i]]<-data.table::fread(file_impute2_info, showProgress = TRUE)
#https://stackoverflow.com/questions/16566799/change-variable-name-in-for-loop-using-r
df_all[[i]]$info<-as.numeric(as.character(df_all[[i]]$info))
###get length of a0 allele 
a0_length<-unlist(lapply(df_all[[i]]$a0, function(x) nchar(x)))

###throw whatver is structual variants
#df_all[[i]]<-df_all[[i]][which(a0_length==1),]

###get length of a1 allele 
a1_length<-unlist(lapply(df_all[[i]]$a1, function(x) nchar(x)))

###throw whatver was an structual variants
df_all[[i]]<-df_all[[i]][which(a1_length==1 & a0_length==1),]

df_all[[i]]$info<-as.numeric(as.character(df_all[[i]]$info))

#--work with polymorphic SNPs
df_all[[i]]$exp_freq_a1<-as.numeric(as.character(df_all[[i]]$exp_freq_a1))
#df_all[[i]]<-df_all[[i]][which(df_all[[i]]$exp_freq_a1!=0),]

df_all[[i]]<-df_all[[i]][,-c(9,10,11,12)] #delete type, certainty and other useless

###Fix alelel frequency
###
indices<-df_all[[i]]$exp_freq_a1>0.5
df_all[[i]]$exp_freq_a1[indices]<-(1-df_all[[i]]$exp_freq_a1)[indices]
}

df_merged2 <- data.table::rbindlist(df_all, idcol = TRUE)
df_merged2$panel<-"1000G"
print("we are done reading info files from 1000G")
print(dim(df_merged2))
#
#Got to HRC data now
#
setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/HRC_WellcomeTrust/imputation_1000indiv/impute2/HRC/merged_data/merged_info_files/")
filename<-"merged_info_CHR"
ext<-".info"
df_all <- as.list(rep("", 22)) 

for(i in seq(1,22)){
file_impute2_info<-paste(filename,i,ext,sep="")

newdf<-paste("df",i,sep="")

df_all[[i]]<-data.table::fread(file_impute2_info, showProgress = TRUE)
#https://stackoverflow.com/questions/16566799/change-variable-name-in-for-loop-using-r

###get length of a0 allele 
a0_length<-unlist(lapply(df_all[[i]]$a0, function(x) nchar(x)))

###throw whatver is structual variants
#df_all[[i]]<-df_all[[i]][which(a0_length==1),]

###get length of a1 allele 
a1_length<-unlist(lapply(df_all[[i]]$a1, function(x) nchar(x)))

###throw whatver was an structual variants
#df_all[[i]]<-df_all[[i]][which(a1_length==1),]
df_all[[i]]<-df_all[[i]][which(a1_length==1 & a0_length==1),]

#--work with polymorphic SNPs
df_all[[i]]$info<-as.numeric(as.character(df_all[[i]]$info))

df_all[[i]]$exp_freq_a1<-as.numeric(as.character(df_all[[i]]$exp_freq_a1))

df_all[[i]]<-df_all[[i]][,-c(9,10,11,12)] #delete type, certainty and other useless

###Fix alelel frequency
###
indices<-df_all[[i]]$exp_freq_a1>0.5
df_all[[i]]$exp_freq_a1[indices]<-(1-df_all[[i]]$exp_freq_a1)[indices]
}

df_merged <- data.table::rbindlist(df_all, idcol = TRUE)
df_merged$panel<-"HRC"

print("we are done reading info files from HRC")
print(dim(df_merged))

####################################
#Merge data frames
####################################

setwd("/mnt/mfs/scratch/GT_Admix/ss5505/imputation/batch2/plots/impute2/")
merged_dfs<-rbind(df_merged2,df_merged)
merged_dfs<-merged_dfs[,-c(1,2,3,4,5,6,9)]
print("done with merging")
print("Object size is ")
print(object.size(merged_dfs),units = "Mb")

binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1, include.lowest = TRUE,seq(0, 1, by=0.05)))
binned_snps<-binned_snps[complete.cases(binned_snps),]

median_all<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_all,"vio_all",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

theme_set(theme_grey(base_size = 18))
CairoJPEG(filename = "vioplot_HRC1000G2.jpeg",quality = 75,width = 1500, height = 900)
ggplot(binned_snps,aes(x=binned,y=info,fill=panel)) + scale_y_continuous(breaks=seq(0.0, 1, 0.2)) +
geom_violin(trim = FALSE) + scale_fill_manual(values=c("grey55","cyan")) +
ggtitle("Impute2: HRC and 1000G") +
xlab("Minor Allele Frequency") + ylab("Info quality") + theme(legend.position="top" , plot.title = element_text(hjust = 0.5)) +
theme(axis.text=element_text(size=17)) + coord_flip()
dev.off()

print("done with All MAF")
binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1, include.lowest = TRUE,seq(0, 0.05, by=0.005)))
binned_snps<-binned_snps[complete.cases(binned_snps),]

median_uncommon<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_uncommon,"vio_median_uncommon",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

CairoJPEG(filename = "vioplot_HRC1000G_uncommon2.jpeg",quality = 75,width = 1500, height = 900)
theme_set(theme_grey(base_size = 18))
ggplot(binned_snps,aes(x=binned,y=info,fill=panel)) +  scale_y_continuous(breaks=seq(0.0, 1, 0.2)) +
geom_violin(trim = FALSE) + scale_fill_manual(values=c("grey55","cyan")) +
ggtitle("Impute2: HRC and 1000G quality within MAF<=5%") +
xlab("Minor Allele Frequency") + ylab("Info quality") + theme(legend.position="top" , plot.title = element_text(hjust = 0.5)) +
theme(axis.text=element_text(size=17)) + coord_flip()
dev.off()

print("done with 5% MAF")

binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1, include.lowest = TRUE,seq(0, 0.01, by=0.001)))
binned_snps<-binned_snps[complete.cases(binned_snps),]

median_rare<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_rare,"vio_median_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

CairoJPEG(filename = "vioplot_HRC1000G_ultrarare2.jpeg",quality = 75,width = 1500, height = 900)
theme_set(theme_grey(base_size = 18))
ggplot(binned_snps,aes(x=binned,y=info,fill=panel)) +  scale_y_continuous(breaks=seq(0.0, 1, 0.2)) +
geom_violin(trim = FALSE) + scale_fill_manual(values=c("grey55","cyan")) +
ggtitle("Impute2: HRC and 1000G quality within MAF<=1%") +
xlab("Minor Allele Frequency") + ylab("Info quality") + theme(legend.position="top" , plot.title = element_text(hjust = 0.5)) +
theme(axis.text=element_text(size=17)) + coord_flip()
dev.off()

print("Done with ultra rare")
