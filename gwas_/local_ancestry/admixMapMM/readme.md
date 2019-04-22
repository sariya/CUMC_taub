#Date 04/22/2019
#
#per chromosome perform admix mapping mixed model

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
[1] parallel  stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
[1] gdsfmt_1.14.1       GWASTools_1.24.1    Biobase_2.38.0
[4] BiocGenerics_0.24.0

loaded via a namespace (and not attached):
 [1] zoo_1.8-5          tidyselect_0.2.5   purrr_0.3.2        reshape2_1.4.3
 [5] DNAcopy_1.52.0     splines_3.4.2      lattice_0.20-38    GWASExactHW_1.01
 [9] mgcv_1.8-28        pan_1.6            blob_1.1.1         survival_2.44-1.1
[13] rlang_0.3.3        jomo_2.6-7         pillar_1.3.1       nloptr_1.2.1
[17] foreign_0.8-71     glue_1.3.1         DBI_1.0.0          bit64_0.9-7
[21] plyr_1.8.4         stringr_1.4.0      MatrixModels_0.4-1 psych_1.8.12
[25] memoise_1.1.0      SparseM_1.77       lmtest_0.9-36      quantreg_5.38
[29] broom_0.4.4        Rcpp_1.0.1         quantsmooth_1.44.0 bit_1.1-14
[33] lme4_1.1-17        mnormt_1.5-5       digest_0.6.18      stringi_1.4.3
[37] dplyr_0.8.0.1      grid_3.4.2         tools_3.4.2        sandwich_2.5-0
[41] magrittr_1.5       tibble_2.1.1       RSQLite_2.1.1      mice_3.4.0
[45] crayon_1.3.4       tidyr_0.8.3        pkgconfig_2.0.2    MASS_7.3-51.3
[49] Matrix_1.2-17      assertthat_0.2.1   minqa_1.2.4        logistf_1.22
[53] rpart_4.1-13       mitml_0.3-7        R6_2.4.0           nnet_7.3-12
[57] nlme_3.1-137       compiler_3.4.2
