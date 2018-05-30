#!/bin/Rscript

###########################################
###Date 05/04/2018
####Sanjeev Sariya 
#Running under: Debian GNU/Linux 9 (stretch) R version 3.4.2
# GWASTools_1.24.1    Biobase_2.38.0      BiocGenerics_0.24.0
# GENESIS_2.10.0      SNPRelate_1.14.0    gdsfmt_1.14.1
##########################################

library("SNPRelate",lib.loc="/mnt/mfs/hgrcgrid/homes/ss5505/R_LIB")
library("GENESIS",lib.loc="/mnt/mfs/hgrcgrid/homes/ss5505/R_LIB")
library("GWASTools",lib.loc="/mnt/mfs/hgrcgrid/homes/ss5505/R_LIB")

bed<-"ordered_merged_common_LD0.3.bed"
fam<-"ordered_merged_common_LD0.3.fam"
bim<-"ordered_merged_common_LD0.3.bim"

snpgdsBED2GDS(bed, fam, bim, "10758.gds")
snpgdsSummary("10758.gds")
print(snpgdsSummary("10758.gds"))

# read in GDS data
geno <- GdsGenotypeReader(filename = "10758.gds")
# create a GenotypeData class object
genoData <- GenotypeData(geno)

iids <- getScanID(genoData)
print(head(iids))

##########
#####make king 
####   king -b 1066_PR_updated_fids_iids_mat_pat_ids.bed --kinship --prefix kinshipKing

KINGmat <- king2mat(file.kin0 = "kinshipKing.kin0", iids = iids, file.kin="kinshipKing.kin")
#
mypcair <- pcair(genoData = genoData, kinMat = KINGmat, divMat = KINGmat)

#--print unrelated
write.table(mypcair$unrels,"unrelated_output",quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

print("we have printed unrelated")

#--print related
write.table(mypcair$rels,"related_output",quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

print("we have printed related")

#--print first ten PCs for entire data
write.table((mypcair$vectors)[,1:10],"first_tenPCs",quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

print("we have printed PCs")

mypcrelate <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:2], training.set = mypcair$unrels, write.to.gds = TRUE)

kinship_dataframe<-mypcrelate$kinship
write.table(kinship_dataframe,"_output_kinship",quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

close(genoData)

#----print other related information about data ----

print("Related people we have")
print(length(mypcair$rels))

print("unRelated people we have")
print(length(mypcair$unrels))

print("unRelated people we have")
print(length(mypcair$unrels))

print("Total people we have")
print(mypcair$nsamp)

print("Total SNPs we have")
print(mypcair$nsnps)
print("We are exiting")