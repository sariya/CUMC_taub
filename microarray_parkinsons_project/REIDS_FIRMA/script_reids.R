#!/bin/Rscript

#
#
#Sanjeev Sariya
#11/13/2018
#Gius - KM PD expression data
#
#
jj<-DataProcessing(chipType = "HuEx-1_0-st-v2",
ExonSummarization = TRUE, 
tags = "coreR3,A20071112,EP",GeneSummarization = TRUE, 
FIRMA = TRUE,
Location = "./",
Name = "CEL", verbose = TRUE)


"
R version 3.4.3 (2017-11-30)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 17134)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] REIDS_0.1.0            GenomeGraphs_1.38.0    biomaRt_2.34.2         aroma.light_3.8.0      aroma.affymetrix_3.1.1 affxparser_1.50.0      aroma.core_3.1.3       R.devices_2.16.0       R.filesets_2.12.1      R.utils_2.7.0         
[11] R.oo_1.22.0            R.methodsS3_1.7.1     

loaded via a namespace (and not attached):
 [1] progress_1.2.0       R.huge_0.9.0         zoo_1.8-4            DNAcopy_1.52.0       listenv_0.7.0        lattice_0.20-35      stats4_3.4.3         base64enc_0.1-3      MCMCpack_1.4-4       blob_1.1.1           XML_3.98-1.16       
[12] rlang_0.3.0.1        DBI_1.0.0            RColorBrewer_1.1-2   BiocGenerics_0.24.0  bit64_0.9-7          matrixStats_0.54.0   R.cache_0.13.0       stringr_1.3.1        MatrixModels_0.4-1   future_1.10.0        codetools_0.2-15    
[23] coda_0.19-2          memoise_1.1.0        Biobase_2.38.0       SparseM_1.77         IRanges_2.12.0       lmtest_0.9-36        Cairo_1.5-9          quantreg_5.36        parallel_3.4.3       AnnotationDbi_1.40.0 Rcpp_0.12.19        
[34] PSCBS_0.64.0         S4Vectors_0.16.0     R.rsp_0.43.0         bit_1.1-14           mcmc_0.9-5           hms_0.4.2            digest_0.6.18        stringi_1.1.7        tools_3.4.3          bitops_1.0-6         magrittr_1.5        
[45] RCurl_1.95-4.11      RSQLite_2.1.1        crayon_1.3.4         future.apply_1.0.1   pkgconfig_2.0.2      Matrix_1.2-12        MASS_7.3-47          data.table_1.11.8    prettyunits_1.0.2    assertthat_0.2.0     httr_1.3.1          
[56] R6_2.3.0             globals_0.12.4       aroma.apd_0.6.0      compiler_3.4.3      

"