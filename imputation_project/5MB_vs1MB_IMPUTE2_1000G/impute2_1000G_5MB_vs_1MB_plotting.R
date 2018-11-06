
#!/bin/Rscript
#--
#Compare 5MB with 1MB imputation 1oooG IMPUTE2
library(plyr)
library(data.table)
library(Cairo)
library(ggplot2)

df_info_impute_1MB<-data.table::fread("1MB_CHR14.info", showProgress = TRUE) # 147923     12
df_info_impute_1MB<-df_info_impute_1MB[,-c(8,9,10,11,12)] #delete type, certainty and other useless

###get length of a0 allele 
a0_length<-unlist(lapply(df_info_impute_1MB$a0, function(x) nchar(x)))

###get length of a1 allele 
a1_length<-unlist(lapply(df_info_impute_1MB$a1, function(x) nchar(x)))

###throw whatver was an structual variants
df_info_impute_1MB<-df_info_impute_1MB[which(a1_length==1 & a0_length ==1),]

#-just check if dim equals sun of length of $a0
dim(df_info_impute_1MB)
sum(nchar(df_info_impute_1MB$a1))

indices<-df_info_impute_1MB$exp_freq_a1>0.5
df_info_impute_1MB$exp_freq_a1[indices]<-(1-df_info_impute_1MB$exp_freq_a1)[indices]


df_info_impute_5MB<-data.table::fread( "5mbCHR14_info.info" , showProgress = TRUE) # 147923     12
df_info_impute_5MB<-df_info_impute_5MB[,-c(8,9,10,11,12)] #delete type, certainty and other useless

###get length of a0 allele 
a0_length<-unlist(lapply(df_info_impute_5MB$a0, function(x) nchar(x)))


###get length of a1 allele 
a1_length<-unlist(lapply(df_info_impute_5MB$a1, function(x) nchar(x)))

###throw whatver was an structual variants
df_info_impute_5MB<-df_info_impute_5MB[which(a1_length==1 & a0_length ==1),]

#-just check if dim equals sun of length of $a0
dim(df_info_impute_5MB)
sum(nchar(df_info_impute_5MB$a1))

indices<-df_info_impute_5MB$exp_freq_a1>0.5
df_info_impute_5MB$exp_freq_a1[indices]<-(1-df_info_impute_5MB$exp_freq_a1)[indices]

df_info_impute_5MB$chunk="5MB"
df_info_impute_1MB$chunk="1MB"
merged_dfs<-rbind(df_info_impute_5MB,df_info_impute_1MB)

#--set MAF range to get average
maf_range<-c(0,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5)
df_maf_range <- as.list(rep("", length(maf_range))) 
x=1
for ( i in maf_range) {

temp_df<-as.data.frame(merged_dfs %>% filter(exp_freq_a1==i) %>% group_by(chunk) %>% summarize_at(vars(info),mean))

if(nrow(temp_df) !=0){
#at times we don't see output for a particular MAF
df_maf_range[[x]]<-temp_df %>%  plyr::rename(c('info'='ave')) %>% mutate(freq=i)

}
x=x+1
}

temp_merge<-rbindlist(lapply(df_maf_range, as.data.table),fill=TRUE)
temp_merge<-temp_merge[,-c(4)]
temp_merge<-temp_merge[complete.cases(temp_merge),]
colnames(temp_merge)<-c("Chunk_length","ave","freq")
write.table(temp_merge,"df_maf_5_1MB_output",col.names=TRUE,quote=FALSE,row.names=FALSE,sep="\t")

##########################
###Plotting
#########################
CairoPDF(file = "5MB_1MBlinegraph.pdf")

df_freq_infoave<-read.table("df_maf_5_1MB_output",header=TRUE)

ggplot(df_freq_infoave,aes(x=as.character(freq),y=ave,color=Chunk_length,group=Chunk_length,linetype=Chunk_length)) + 
scale_y_continuous(breaks=seq(0.0, 1, 0.2)) +geom_line(size=1) + 
geom_point(size = 1.5) + ggtitle("Comparison of average info for 5MB vs 1MB") +
xlab("Minor Allele Frequency") + ylab("Average info quality") + theme(legend.position="right" , plot.title = element_text(hjust = 0.5)) +
theme(axis.text=element_text(size=9,family="Times New Roman")) +  theme(plot.title = element_text(size=13,family="Times New Roman"))

dev.off()


#jpeg("linegrapxx5mbvs1mb_1000G.jpeg",quality = 75,width = 1500, height = 900)
#ggplot(df_freq_infoave,aes(x=freq,y=ave,color=Panel,group=Panel,linetype=Panel)) + 
#scale_y_continuous(breaks=seq(0.0, 1, 0.2)) +geom_line(size=1) + 
#geom_point(size = 1.5) + ggtitle("Comparison for average info between HRC and 1000G") +
#xlab("Minor Allele Frequency") + ylab("Average info quality") + theme(legend.position="right" , plot.title = element_text(hjust = 0.5)) +
#theme(axis.text=element_text(size=11)) +  theme(plot.title = element_text(size=14))
#dev.off()
#dev.off()


#jpeg("xx5mbvs1mb_1000G.jpeg",quality = 75,width = 1500, height = 900)
#ggplot(df_info_impute_5MB, aes(x = "5MB" , y = info)) + geom_boxplot (alpha=0.2,notch=FALSE,width=0.1,fill="#56B4E9", outlier.colour = "red",outlier.size = 1) +
#geom_boxplot (data=df_info_impute_1MB,aes(x = "1MB" , y = info),alpha=0.2,notch=FALSE,width=0.1,fill="Red", outlier.colour = "red",outlier.size = 1) +
#ggtitle("CHR14: 5MB vs 1MB: IMPUTE2 1000G") + xlab("Tool") + ylab("Info") + 
#theme(legend.position="top" , plot.title = element_text(hjust = 0.4)) +
#theme(axis.text=element_text(size=18))
#dev.off()
#

#
#Check for difference in quality in rare variants
#

#df_info_impute_5MB$chunk="5MB"
#df_info_impute_1MB$chunk="1MB"
#merged_dfs<-rbind(df_info_impute_5MB,df_info_impute_1MB)

#
## binned_snps<-transform(df_info_impute_1MB,binned=cut(df_info_impute_1MB$exp_freq_a1, include.lowest = TRUE,seq(0, 1, by=0.05)))
# jj<-as.data.frame(transform(df_info_impute_1MB,binned=cut(df_info_impute_1MB$exp_freq_a1, include.lowest = TRUE,seq(0, 1, by=0.05))) %>% group_by(binned,chunk) %>% summarize_at(vars(info),median))
#mmm <- melt(jj)
#ggplot(data = mmm, aes(x =  binned , y = value, fill = chunk)) +        geom_bar(stat = 'identity', position = 'dodge')

#binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 1, by=0.10)))
#binned_snps<-binned_snps[complete.cases(binned_snps),]

#
#All variants
#
#jpeg("5mbvs1mb_1000G_Binned_all.jpeg",quality = 75,width = 1500, height = 900)
#theme_set(theme_grey(base_size = 18))
#ggplot(binned_snps,aes(x=binned,y=info,fill=chunk)) + scale_y_continuous(breaks=seq(0.0, 1, 0.1)) +
#geom_boxplot( alpha=0.2,notch=FALSE,width=0.5) + scale_fill_manual(values=c("blue","green")) +
#ggtitle("Impute2: 5MB vs 1MB All variants") +
#xlab("Minor Allele Frequency") + ylab("Info quality") + theme(legend.position="top" , plot.title = element_text(hjust = 0.5)) +
#theme(axis.text=element_text(size=17)) + coord_flip()
#dev.off()

#
#5%MAF
#

#binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 0.05, by=0.005)))
#binned_snps<-binned_snps[complete.cases(binned_snps),]

#jpeg("5mbvs1mb_1000G_maf5Perc.jpeg",quality = 75,width = 1500, height = 900)
#theme_set(theme_grey(base_size = 18))
#ggplot(binned_snps,aes(x=binned,y=info,fill=chunk)) + scale_y_continuous(breaks=seq(0.0, 1, 0.1)) +
#geom_boxplot( alpha=0.2,notch=FALSE,width=0.5) + scale_fill_manual(values=c("blue","green")) +
#ggtitle("Impute2: 5MB vs 1MB MAF<=5%") +
#xlab("Minor Allele Frequency") + ylab("Info quality") + theme(legend.position="top" , plot.title = element_text(hjust = 0.5)) +
#theme(axis.text=element_text(size=17)) + coord_flip()
#dev.off()

#
#Ultra rare
#

#binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 0.01, by=0.001)))
#binned_snps<-binned_snps[complete.cases(binned_snps),]


#jpeg("5mbvs1mb_1000G_ultrarare.jpeg",quality = 75,width = 1500, height = 900)
#theme_set(theme_grey(base_size = 18))
#ggplot(binned_snps,aes(x=binned,y=info,fill=chunk)) + scale_y_continuous(breaks=seq(0.0, 1, 0.1)) +
#geom_boxplot( alpha=0.2,notch=FALSE,width=0.5) + scale_fill_manual(values=c("blue","green")) +
#ggtitle("Impute2: 5MB vs 1MB MAF<=1%") +
#xlab("Minor Allele Frequency") + ylab("Info quality") + theme(legend.position="top" , plot.title = element_text(hjust = 0.5)) +
#theme(axis.text=element_text(size=17)) + coord_flip()
#dev.off()


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
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] ggplot2_3.0.0     Cairo_1.5-9       data.table_1.11.8 plyr_1.8.4

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.19     withr_2.1.2      assertthat_0.2.0 crayon_1.3.4
 [5] dplyr_0.7.7      grid_3.4.2       R6_2.3.0         gtable_0.2.0
 [9] magrittr_1.5     scales_1.0.0     pillar_1.3.0     rlang_0.3.0
[13] lazyeval_0.2.1   bindrcpp_0.2.2   glue_1.3.0       purrr_0.2.5
[17] munsell_0.5.0    compiler_3.4.2   pkgconfig_2.0.2  colorspace_1.3-2
[21] tidyselect_0.2.5 bindr_0.1.1      tibble_1.4.2

"

