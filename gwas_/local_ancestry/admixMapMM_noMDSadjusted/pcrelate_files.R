#!/usr/bin/Rscript

#
#Date 04/23/2019
#Sanjeev Sariya
#Local Ancestry
#https://rdrr.io/bioc/GENESIS/man/admixMapMM.html
#
#We need pcrelate data for admixmapMM. Use this script for .gds file creation in intermediatery steps
#
print("Begin Script")

library("SNPRelate")
library("GENESIS")
library("GWASTools")

print("loaded libraries")

bed<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/genesis_admixMapMM/sorted_merged_updatedIDs_commonSNPs.bed"
fam<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/genesis_admixMapMM/sorted_merged_updatedIDs_commonSNPs.fam"
bim<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/genesis_admixMapMM/sorted_merged_updatedIDs_commonSNPs.bim"


snpgdsBED2GDS(bed, fam, bim, "10758.gds")
snpgdsSummary("10758.gds")
print(snpgdsSummary("10758.gds"))

# read in GDS data
geno <- GdsGenotypeReader(filename = "10758.gds")
# create a GenotypeData class object
genoData <- GenotypeData(geno)

iids <- getScanID(genoData)
print(head(iids))

write.table(iids,"sample.id",quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

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

########################################use output in admixmap from following step
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