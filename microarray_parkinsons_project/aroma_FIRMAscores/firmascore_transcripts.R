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

###
###Prove level annotation
###

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
affy <- "affy_huex_1_0_st_v2"
genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "affy_huex_1_0_st_v2"), filters = "affy_huex_1_0_st_v2", values=c("2315101","2315102","2315103","2315104"), mart = mart)

