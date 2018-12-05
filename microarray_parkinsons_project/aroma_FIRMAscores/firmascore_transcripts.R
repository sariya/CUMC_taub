#
#Date 04 11 2018
# Sanjeev Sariya transcript level abundances

#
#http://www.aroma-project.org/vignettes/FIRMA-HumanExonArrayAnalysis/
##All analyses are performed on Windows

library("aroma.affymetrix")

cdf <- AffymetrixCdfFile$byChipType("HuEx-1_0-st-v2", tags="coreR3,A20071112,EP")
print(cdf)
cs <- AffymetrixCelSet$byName("alldata", cdf=cdf) # rawdata in cwd. HuEx-1_0-st-v2 in celfiles # this should match cdf file names too
print(cs)
setCdf(cs, cdf)
bc <- RmaBackgroundCorrection(cs) #backgroun
csBC <- process(bc,verbose=verbose)

qn <- QuantileNormalization(csBC, typesToUpdate="pm") 
csN <- process(qn) 

getCdf(csN) ##get CDF infor
setCdf(csN, cdf)
plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE,tag=c("*","coreR3")) #entrire transcript 
print(plmTr)

fit(plmTr, verbose=verbose,force=F)

#----alternative splicing for transcript PLM
firma <- FirmaModel(plmTr,tag=c("*","coreR3")) #The FIRMA analysis only works from the PLM based on transcripts.
fit(firma, verbose=verbose)
fs <- getFirmaScores(firma)
cestrfs <- extractDataFrame(fs, units=NULL, addNames=TRUE) #extract the FIRMA scores -- but Nan
head(cestrfs)
print(dim(cestrfs))

#--write as is
write.table(cestrfs,"firmascores",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

temp_firmacolnames<-(colnames(cestrfs))[6:length(colnames(cestrfs))]

#--fix col names

for(i in 1:length(temp_firmacolnames)){

temp_firmacolnames[i]<-(paste("X",temp_firmacolnames[i],sep=""))
temp_firmacolnames[i]<-(paste(temp_firmacolnames[i],".CEL",sep=""))
temp_firmacolnames[i]<-gsub("-",".",temp_firmacolnames[i])
}

temp_firmacolnames<-c(colnames(cestrfs)[1:5],temp_firmacolnames)
print(length(temp_firmacolnames))

colnames(cestrfs)<-temp_firmacolnames
#--write as column fixed

write.table(cestrfs,"firmascores_colnamesfixed",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

list_split_firmaoutput<-split(cestrfs, (seq(nrow(cestrfs))-1) %/% 30000) 
print(length(list_split_firmaoutput))

#
#Split into multiple chunks
#

for(i in 1:length(list_split_firmaoutput) )
{
file<-paste("aroma.firmascore_chunk",i,sep="_") #--make file name for chunk
file_dest<-paste("./",file,sep="/")
write.table(list_split_firmaoutput[[i]], file_dest, sep="\t", quote=F, row.names=FALSE,col.names=TRUE)
}

###
###Prove level annotation
###

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
affy <- "affy_huex_1_0_st_v2"
genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "affy_huex_1_0_st_v2"), filters = "affy_huex_1_0_st_v2", values=c("2315101","2315102","2315103","2315104"), mart = mart)

"
R version 3.4.3 (2017-11-30)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 17134)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] preprocessCore_1.40.0  aroma.light_3.8.0      aroma.affymetrix_3.1.1 affxparser_1.50.0      aroma.core_3.1.3       R.devices_2.16.0       R.filesets_2.12.1      R.utils_2.7.0          R.oo_1.22.0            R.methodsS3_1.7.1     

loaded via a namespace (and not attached):
 [1] DNAcopy_1.52.0       zlibbioc_1.24.0      BiocGenerics_0.24.0  aroma.apd_0.6.0      R.cache_0.13.0       globals_0.12.4       tools_3.4.3          parallel_3.4.3       Biobase_2.38.0       R.huge_0.9.0         affy_1.56.0         
[12] PSCBS_0.64.0         matrixStats_0.54.0   digest_0.6.18        affyio_1.48.0        R.rsp_0.43.0         base64enc_0.1-3      codetools_0.2-15     future.apply_1.0.1   BiocInstaller_1.28.0 compiler_3.4.3       future_1.10.0       
[23] listenv_0.7.0        Cairo_1.5-9         


"