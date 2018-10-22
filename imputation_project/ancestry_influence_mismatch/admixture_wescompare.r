
#!/bin/Rscript

#Date 07/19/2018
#

#
#Take mismatch from WES 14th chromosomes from VCF-Compare. 
#
#Out of 1,000 individuals 262 are WESed
#Find common SNPs on 14th Chromosome WES-HRC-1000G. These SNPs must have 0% missingness in WES data
#Convert genotype probabilities from IMPUTE2 of these common SNPS into plink files.
#Make VCF files, compare on different info threshold and MAF bins. Get mismatch counts for Info80%, MAF1-5%, MAF0.1-1% and MAF 0-0.1%
# The mismtach are from WES-HRC comparison and HRC-WES for SNPs that are 100% present in WES data, overlapping in HRC-1000G
#Get Admixture component for 262 people for which WES-HRC-1000G imputation are done. 
#
#Run poisson regression using mismatches as outcome and Global ancestry as predictors
##


library(dplyr)

#global ancestry file
df_global<-data.table::fread("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/accuracy_imputation/gius_accuracy/info_maf_bins/global_admixture", showProgress = TRUE,header=TRUE)

#wes file mistmatch per person
file.mismatch<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/accuracy_imputation/gius_accuracy/info_maf_bins/10.19.2018_analysis/HRC_100perc/info0.8/MAF_0_0.1PERC/vcf_compare_0_0.1perc"

df_ref_wes<-data.table::fread(file.mismatch, showProgress = TRUE)

#do a left join
join_ref<-left_join(df_global,df_ref_wes,by=c("IID"="V1"))

#remove incomplete cases
join_ref<-join_ref[complete.cases(join_ref),]

print( dim(join_ref))

summary(glm(V2 ~ CEU, family="poisson", data=join_ref))$coefficients[,4]
#--slope
coef((glm(V2 ~ CEU, family="poisson", data=join_ref)))["CEU"]
#--intercept
coef((glm(V2 ~ CEU, family="poisson", data=join_ref)))["(Intercept)"]
#--Pvalue
summary(glm(V2 ~ CEU, family="poisson", data=join_ref))$coeff[2,4]

summary(glm(V2 ~ NAT, family="poisson", data=join_ref))$coefficients[,4]

#--slope
coef((glm(V2 ~ NAT, family="poisson", data=join_ref)))["NAT"]
#--intercept
coef((glm(V2 ~ NAT, family="poisson", data=join_ref)))["(Intercept)"]
#--Pvalue
summary(glm(V2 ~ NAT, family="poisson", data=join_ref))$coeff[2,4]

summary(glm(V2 ~ YRI, family="poisson", data=join_ref))$coefficients[,4]
#--slope
coef((glm(V2 ~ YRI, family="poisson", data=join_ref)))["YRI"]
#--intercept
coef((glm(V2 ~ YRI, family="poisson", data=join_ref)))["(Intercept)"]
#--Pvalue
summary(glm(V2 ~ YRI, family="poisson", data=join_ref))$coeff[2,4]

"
 sessionInfo()
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
[1] dplyr_0.7.7

loaded via a namespace (and not attached):
 [1] tidyselect_0.2.5 compiler_3.4.2   magrittr_1.5     assertthat_0.2.0
 [5] R6_2.3.0         pillar_1.3.0     bindrcpp_0.2.2   glue_1.3.0
 [9] tibble_1.4.2     crayon_1.3.4     Rcpp_0.12.19     pkgconfig_2.0.2
[13] rlang_0.2.2      purrr_0.2.5      bindr_0.1.1
>


"