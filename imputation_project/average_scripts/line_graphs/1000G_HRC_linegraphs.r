#!/bin/Rscript

##
#Date 06/26/2018
#
#MAF line graph for KGP and HRC
#Find counts per different info thresholds with different MAF bins

#Read KGP data first. 
#

library(dplyr)
library("data.table")
library("Cairo")
library(plyr)
library(ggplot2)

setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/impute2/1KGP/merged_files/unharmonized/info/")
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
df_all[[i]]<-df_all[[i]][which(a1_length==1 & a0_length==1),]
df_all[[i]]$info<-as.numeric(as.character(df_all[[i]]$info))

#--work with polymorphic SNPs
df_all[[i]]$exp_freq_a1<-as.numeric(as.character(df_all[[i]]$exp_freq_a1))

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

setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/impute2/HRC/merged_data/merged_info_files/")

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

###get length of a1 allele 
a1_length<-unlist(lapply(df_all[[i]]$a1, function(x) nchar(x)))

###throw whatver was an structual variants
df_all[[i]]<-df_all[[i]][which(a1_length==1 & a0_length==1),]
df_all[[i]]$info<-as.numeric(as.character(df_all[[i]]$info))

#--work with polymorphic SNPs
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
merged_dfs<-rbind(df_merged2,df_merged)
merged_dfs<-merged_dfs[,-c(1,2,3,4,5,6,9)]
print("done with merging")
print(head(merged_dfs))
print(tail(merged_dfs))
setwd("/mnt/mfs/scratch/GT_Admix/ss5505/imputation/batch2/plots/impute2/line_graphs/")

#--define range of MAF to calculate Info average

maf_range<-c(0,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5)
df_maf_range <- as.list(rep("", length(maf_range))) 
x=1
for ( i in maf_range) {

temp_df<-as.data.frame(merged_dfs %>% filter(exp_freq_a1==i) %>% group_by(panel) %>% summarize_at(vars(info),mean))

if(nrow(temp_df) !=0){
#at times we don't see output for a particular MAF
df_maf_range[[x]]<-temp_df %>%  rename(c('info'='ave')) %>% mutate(freq=i)

}

##df_maf_range[[x]]<- as.data.frame(merged_dfs %>% filter(exp_freq_a1==i) %>% group_by(panel) %>% summarize_at(vars(info),mean)) %>%  rename(c('info'='ave')) %>% mutate(freq=i)

x=x+1
}

#--there are emplty df in list
#https://stackoverflow.com/a/16591717/2740831
#
temp_merge<-rbindlist(lapply(df_maf_range, as.data.table),fill=TRUE)
temp_merge<-temp_merge[,-c(4)]
temp_merge<-temp_merge[complete.cases(temp_merge),]
#df_freq_infoave <- data.table::rbindlist(df_maf_range)
df_freq_infoave<-temp_merge
df_freq_infoave$freq<-as.character(df_freq_infoave$freq)

write.table(df_freq_infoave,"info_maf_output",col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

CairoJPEG(filename = "linegraph.jpeg",quality = 75,width = 1500, height = 900)
theme_set(theme_grey(base_size = 18))
ggplot(df_freq_infoave,aes(x=freq,y=ave,color=panel,group=panel,linetype=panel)) + 
scale_y_continuous(breaks=seq(0.0, 1, 0.2)) +geom_line(size=1) + 
geom_point(size = 1.5) + ggtitle("Impute2: HRC and 1000G Average Info") +
xlab("Minor Allele Frequency") + ylab("Average Info quality") + theme(legend.position="top" , plot.title = element_text(hjust = 0.5)) +
theme(axis.text=element_text(size=17))
dev.off()


print("Check output or current working directory")

#
#work with bins and calculate mean info and plot 
#
#

uncommon<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.05),]
rare<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.01),]
ultra_rare<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.001),]

write.table(as.data.frame(uncommon %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninfo_uncommon",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(rare %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninfo_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(ultra_rare %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninfo_ultra_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

uncommon_info0.40<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.05 &  merged_dfs$info>=0.40 ),]
rare_info0.40<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.01   & merged_dfs$info>=0.40 ),]
ultra_rare_info0.40<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.001  & merged_dfs$info>=0.40 ),]

write.table(as.data.frame(uncommon_info0.40 %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninfo0.40_uncommon",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(rare_info0.40 %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninfo0.40_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(ultra_rare_info0.40 %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninfo0.40_ultra_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

uncommon_info0.80<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.05 &  merged_dfs$info>=0.80 ),]
rare_info0.80<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.01   & merged_dfs$info>=0.80 ),]
ultra_rare_info0.80<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.001  & merged_dfs$info>=0.80 ),]

write.table(as.data.frame(uncommon_info0.80 %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninfo0.80_uncommon",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(rare_info0.80 %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninfo0.80_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(ultra_rare_info0.80 %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninfo0.80_ultra_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")


########
########
########

create_binned_df<-function(df_input,lower_limit,upper_limit,bin_value){

#--function to create dfs based on lower, upper lomit provided by used. and breaks we need.

return (transform(df_input,binned=cut(df_input$exp_freq_a1,include.lowest = TRUE, seq(lower_limit, upper_limit, by=bin_value))))
}

########
########
########
create_plots<-function(img_name,x_lower,x_upper){
print("plot ")
}

binned_snps<-create_binned_df(merged_dfs,0,1,0.05)
#binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 1, by=0.05)))
binned_snps<-binned_snps[complete.cases(binned_snps),]

mean_all<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),mean))
write.table(mean_all,"mean_all",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

binned_snps<-create_binned_df(merged_dfs,0,0.01,0.001)  ###transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 0.01, by=0.001)))
binned_snps<-binned_snps[complete.cases(binned_snps),]
mean_rare<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),mean))
write.table(mean_rare,"mean_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

#
#Work with MAF ranges
#
#

df1_5MAF<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.05 & merged_dfs$exp_freq_a1>=0.01),]
df0.1_1MAF<-merged_dfs[which(merged_dfs$exp_freq_a1<0.01 & merged_dfs$exp_freq_a1>=0.001 ),]
df0_0.1MAF<-merged_dfs[which(merged_dfs$exp_freq_a1<0.001 & merged_dfs$exp_freq_a1>=0),]

write.table(as.data.frame(df1_5MAF %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninforange_1_5MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(df0.1_1MAF %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninforange_0.1_1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(df0_0.1MAF %>% group_by(panel) %>% summarize_at(vars(info),mean)),"meaninforange_0_0.1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

write.table(as.data.frame(df1_5MAF %>% group_by(panel) %>% summarize_at(vars(info),length)),"countinfo0inforange_1_5MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(df0.1_1MAF %>% group_by(panel) %>% summarize_at(vars(info),length)),"countinfo0inforange_0.1_1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(df0_0.1MAF %>% group_by(panel) %>% summarize_at(vars(info),length)),"countinfo0inforange_0_0.1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")


#####################################################################################################################################
#
#work with Info of 0.40
#####################################################################################################################################

info_1_5MAF<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.05 & merged_dfs$exp_freq_a1>=0.01 & merged_dfs$info>=0.40),]
info_0.1_1MAF<-merged_dfs[which(merged_dfs$exp_freq_a1<0.01 & merged_dfs$exp_freq_a1>=0.001  & merged_dfs$info>=0.40),]
info_0_0.1MAF<-merged_dfs[which(merged_dfs$exp_freq_a1<0.001 & merged_dfs$exp_freq_a1>=0  & merged_dfs$info>=0.40),]

write.table(as.data.frame(info_1_5MAF %>% group_by(panel) %>% summarize_at(vars(info),mean)),"mean0.4inforange_1_5MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_0.1_1MAF %>% group_by(panel) %>% summarize_at(vars(info),mean)),"mean0.4inforange_0.1_1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_0_0.1MAF %>% group_by(panel) %>% summarize_at(vars(info),mean)),"mean0.4inforange_0_0.1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

#
#Get counts per panel
#

write.table(as.data.frame(info_1_5MAF %>% group_by(panel) %>% summarize_at(vars(info),length)),"counts0.4inforange_1_5MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_0.1_1MAF %>% group_by(panel) %>% summarize_at(vars(info),length)),"counts0.4inforange_0.1_1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_0_0.1MAF %>% group_by(panel) %>% summarize_at(vars(info),length)),"countsn0.4inforange_0_0.1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")


#####################################################################################################################################
#
#work with Info of 0.80
#####################################################################################################################################

info_1_5MAF_0.80<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.05 & merged_dfs$exp_freq_a1>=0.01 & merged_dfs$info>=0.80),]
info_0.1_1MAF_0.80<-merged_dfs[which(merged_dfs$exp_freq_a1<0.01 & merged_dfs$exp_freq_a1>=0.001  & merged_dfs$info>=0.80),]
info_0_0.1MAF_0.80<-merged_dfs[which(merged_dfs$exp_freq_a1<0.001 & merged_dfs$exp_freq_a1>=0  & merged_dfs$info>=0.80),]

write.table(as.data.frame(info_1_5MAF_0.80 %>% group_by(panel) %>% summarize_at(vars(info),mean)),"mean0.8inforange_1_5MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_0.1_1MAF_0.80 %>% group_by(panel) %>% summarize_at(vars(info),mean)),"mean0.8inforange_0.1_1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_0_0.1MAF_0.80 %>% group_by(panel) %>% summarize_at(vars(info),mean)),"mean0.8inforange_0_0.1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")


#
#Get counts per panel for 0.80
#

write.table(as.data.frame(info_1_5MAF_0.80 %>% group_by(panel) %>% summarize_at(vars(info),length)),"counts0.8inforange_1_5MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_0.1_1MAF_0.80 %>% group_by(panel) %>% summarize_at(vars(info),length)),"counts0.8inforange_0.1_1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_0_0.1MAF_0.80 %>% group_by(panel) %>% summarize_at(vars(info),length)),"counts0.8inforange_0_0.1MAF",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")


