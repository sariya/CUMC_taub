#!/usr/bin/Rscript

#
#Date 04/29/2019
#Sanjeev Sariya
#Local Ancestry
#Analyse only for females. Model 1

#https://rdrr.io/bioc/GENESIS/man/admixMapMM.html
print("Begin Script")
library(GWASTools)
library(gdsfmt)
library(dplyr)
library(GENESIS)
library(argparse)
library(SNPRelate)
print("loaded libraries")
parser <- ArgumentParser(description="perform admix MM") ##
parser$add_argument('-c',"--ceu",help="CEU input file",required=TRUE) ##CEU Rfmix V1 components
parser$add_argument('-t',"--nat",help="NAT input file",required=TRUE) ## NAT Rfmix V1 components##hgdp
parser$add_argument('-y',"--yri",help="YRI/AFR input file",required=TRUE)  ##YRI RFmix components
parser$add_argument('-x',"--pre",help="prefix for output file",required=TRUE) ##prefix for output
parser$add_argument('-s',"--snps",help="input SNPs file",required=TRUE) # store input snps file
parser$add_argument('-n',"--chr",help="input chr number",required=TRUE) #store chr number

args <- parser$parse_args() #make it a data structure

#
#Set file names
#
#file.afr<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/merged_CHR22_YRI.txt"
#file.nat<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/merged_CHR22_NAT.txt"
#file.ceu<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/merged_CHR22_CEU.txt"
#file.snps<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/CHR22_snps_rfmix"

##make sure sex is M/F. 0/1 isn't acccepted
file.pheno<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/genesis_admixMapMM/input.phenomodel1_onlyfemales"

##kinship is created using GEMMA. Use as is with no other fuss. please. use for all CHRs.
file.kinship<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/kinship/result.cXX.txt"

###external file. use for all CHRs.
file.sampleid<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/genesis_admixMapMM/input.sampleids"

##while workign with RFMix, this file was created per CHR. Use once for all
file.rsPos<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/rsids_allchrs_positions"

##Normalize paths in argparse##
file.afr<-normalizePath(args$yri)
file.nat<-normalizePath(args$nat)
file.ceu<-normalizePath(args$ceu)
file.snps<-normalizePath(args$snps)
chr<-args$chr
out_prefix<-args$pre

print(file.afr)
print(file.nat)
print(file.ceu)
print(file.snps)
print(chr)
print(out_prefix)
#
#Read in files and get their dimensions
#

##chr<-22
chr<-as.integer(as.numeric(as.character(chr))) ##make change to character chromosome number
print(chr)
df.snps<-read.csv(file.snps,header=FALSE)
df.afr<-read.csv(file.afr,header=FALSE)
df.nat<-read.csv(file.nat,header=FALSE)
df.ceu<-read.csv(file.ceu,header=FALSE)

df.sampleids<-read.table(file.sampleid,header=FALSE)
print(dim(df.sampleids))

grm_all <- as.matrix(read.table(file.kinship))
print(dim(grm_all))

rownames(grm_all)<-df.sampleids$V1
colnames(grm_all)<-df.sampleids$V1

df.pheno<-read.table(file.pheno,header=TRUE)
print(dim(df.pheno))
print(dim(df.afr))
print(dim(df.nat))
print(dim(df.ceu))
print(dim(df.snps))

colnames(df.snps)<-c("rsids","a1","a2")
#
#Check nrow in order to ensure data look good.
#		
if(nrow(df.afr)!=nrow(df.nat) | nrow(df.ceu)!=nrow(df.nat)){
    stop("issue with nrow nat ceu yri")
}

if(nrow(df.afr)!=nrow(df.nat) | nrow(df.ceu)!=nrow(df.snps)){
    stop("issue with SNPs nrow")
}

##Check and reading ends#

##Read file with rs ids and chr:pos#

df.rsPos<-read.table(file.rsPos,header=FALSE)
print(dim(df.rsPos))
df.rsPos$V1<-as.character(df.rsPos$V1)
df.rsPos$V2<-as.character(df.rsPos$V2)

temp<-strsplit(df.rsPos$V1,"\\:") #split column for RS
rs_chrmatrix <- do.call(rbind,temp) #do binding
rs_chrmatrix <- as.data.frame(rs_chrmatrix,stringsAsFactors = FALSE)
rs_chrmatrix$rsids<- df.rsPos$V2
rs_chrmatrix$V1<-as.numeric(as.character(rs_chrmatrix$V1))
rs_chrmatrix$V2<-as.numeric(as.character(rs_chrmatrix$V2))

##get length of positions and RSids that match given CHRs#
print(length(which(rs_chrmatrix$V1==chr)))

rs_chrmatrix<-rs_chrmatrix[which(rs_chrmatrix$V1==chr),]

colnames(rs_chrmatrix)<-c("CHR","POS","rsids")
jointed_rschrpos<-left_join(df.snps,rs_chrmatrix,by=c("rsids"))

if(sum(!complete.cases(jointed_rschrpos)) !=0){
    print(chr)

    stop("we have some issue here with Rsids and positions sum not equal to zero!!!!!")
}
##If ends

if(nrow(jointed_rschrpos[complete.cases(jointed_rschrpos),]) != nrow(df.ceu)){
    print(chr)
    stop("we have issue when taking complete cases and nrow with CEU")
}
##if ends

##First three columns are useless in NAT/CEU/AFR #
print(dim(jointed_rschrpos))

print("All data look Good")

print("creating GDS file")

file.gdsdosage<-paste("dosage","_CHR",chr,".gds",sep="")

#snpgdsCreateGeno("test2.gds", 
 #   sample.id = df.sampleids$V1, snp.id = jointed_rschrpos$rsids,
  #  snp.chromosome = jointed_rschrpos$CHR,   # snp.position = jointed_rschrpos$POS)

gdsfile <- createfn.gds(file.gdsdosage)

add.gdsn(gdsfile,"snp.chromosome",rep(chr, nrow(df.ceu)))
#https://github.com/zhengxwen/gdsfmt/issues/23#issuecomment-486452705

##snpgdsCreateGeno(gdsfile ,"sample.id" ,df.sampleids$V1, snp.id = df.ceu[,c(1)],    snp.chromosome = chr,    snp.position = jointed_rschrpos$POS)

matrix.afr<-as.matrix(df.afr[,-c(1:3)]) ##delete first three columns of RS ids and A1 and A2
matrix.nat<-as.matrix(df.nat[,-c(1:3)])
matrix.ceu<-as.matrix(df.ceu[,-c(1:3)])
add.gdsn(gdsfile , "dosage_ceu",matrix.ceu)
add.gdsn(gdsfile , "dosage_nat",matrix.nat)
add.gdsn(gdsfile , "dosage_afr",matrix.afr)

add.gdsn(gdsfile,"sample.id" ,df.sampleids$V1)
add.gdsn(gdsfile,"snp.id" ,df.ceu[,c(1)])
add.gdsn(gdsfile,"snp.position",jointed_rschrpos$POS)

closefn.gds(gdsfile)
##print and get details for verification
gds <- openfn.gds(file.gdsdosage)

id <- read.gdsn(index.gdsn(gds, "snp.id"))
print(length(id ))
chrom <- read.gdsn(index.gdsn(gds, "snp.chromosome"))
print(length(chrom))
closefn.gds(gds) ##close file. :-|

print("Closed and printed information about GDS created")
print("completed creating GDS file")

# fit the null mixed model
covariates <- c( "age") ##make sure sex is M/F. 0/1 isn't acccepted
outcome <- "AD"
# make ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = rownames(grm_all),df.pheno, stringsAsFactors=FALSE))
print("Fit null model")
nullmod <- fitNullMM(scanData = scanAnnot,outcome = outcome,covars = covariates, covMatList = grm_all,verbose = TRUE)
print("Completed fitting null model")

ancestries <- c("nat","afr") ##only two ancestries. 
genoDataList <- list()

tempgds <- openfn.gds(file.gdsdosage) ##open gds dosage and close please
for (anc in  ancestries ){
	print(anc )
  gdsr <- GdsGenotypeReader(tempgds , genotypeVar=paste0("dosage_", anc))
  genoDataList[[anc]] <- GenotypeData(gdsr, scanAnnot=scanAnnot)
}
##don't close the object of GDS, otherwise it gives error Error in ls.gdsn(object@handler) :
##  Invalid GDS node object (it was closed or deleted).

assoc.admix <- admixMapMM(genoDataList,nullMMobj = nullmod)
print(dim(assoc.admix))

##sort based on joint-pval
assoc.admix<- assoc.admix[order(assoc.admix$Joint.pval),] 

##do a merge of P-value output with RSids positions.
assoc.admix.jointed<-left_join(assoc.admix,rs_chrmatrix[,-c(1)],by=c("snpID"="rsids"))

write.table(assoc.admix.jointed, file = out_prefix, append = FALSE, quote = FALSE, sep = "\t", eol = "\n",
            na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape",
        "double"), fileEncoding = "")

closefn.gds(tempgds)

print("Ending Script")
##Important links:
##https://bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#create-a-gwas-snp-gds-file
#
