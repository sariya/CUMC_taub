#!/bin/bash

cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/ratio_exon_to_gene

output=/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/ratio_exon_to_gene/RIN_adjusted/
for i in {1..11}

do

set -x 
file=/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/ratio_exon_to_gene/exon_chunk_ratio_"$i"

#qsub -cwd -N "M1exon$i" -l h_vmem=5G -b y /mnt/mfs/cluster/bin/R-3.4/bin/Rscript model1_analysis.R -e $file -n "$i" -o $output"Model1/" -x "mod1_chunk$i"
#/mnt/mfs/cluster/bin/R-3.4/bin/Rscript model1_analysis.R -e $file -n "$i" -o $output"Model1/" -x "mod1_chunk$i"

#qsub -cwd -N "M2exn$i" -l h_vmem=5G -b y /mnt/mfs/cluster/bin/R-3.4/bin/Rscript model2_analysis.R -e $file -n "$i" -o $output"Model2/" -x "mod2_chunk$i" 
#/mnt/mfs/cluster/bin/R-3.4/bin/Rscript model2_analysis.R -e $file -n "$i" -o $output"Model2/" -x "mod2_chunk$i" 

qsub -cwd -N "M3exn$i" -l h_vmem=5G -b y  /mnt/mfs/cluster/bin/R-3.4/bin/Rscript model3_analysis.R -e $file -n "$i" -o $output"Model3/" -x "mod3_chunk$i"

#qsub -cwd -N "M4exn$i" -l h_vmem=5G -b y  /mnt/mfs/cluster/bin/R-3.4/bin/Rscript model4_compare_cat3_cat4.R -e $file -n "$i" -o $output"compare_cat3_cat4/" -x "mod4_chunk$i"

#qsub -cwd -N "M4exn$i" -l h_vmem=5G -b y /mnt/mfs/cluster/bin/R-3.4/bin/Rscript model4_compare_cat2_cat4.R -e $file -n "$i" -o $output"compare_cat2_cat4/" -x "mod4_chunk$i"
#qsub -cwd -N "M4exn$i" -l h_vmem=5G -b y /mnt/mfs/cluster/bin/R-3.4/bin/Rscript model4_compare_cat2_cat3.R -e $file -n "$i" -o $output"compare_cat2_cat3/" -x "mod4_chunk$i"
#qsub -cwd -N "M4exn$i" -l h_vmem=5G -b y /mnt/mfs/cluster/bin/R-3.4/bin/Rscript model4_compare_cat0_cat4.R -e $file -n "$i" -o $output"compare_cat0_cat4/" -x "mod4_chunk$i"

#qsub -cwd -N "M4exn$i" -l h_vmem=2G -b y /mnt/mfs/cluster/bin/R-3.4/bin/Rscript model4_analysis.R -e $file -n "$i" -o $output -x "mod4_allCat_chunk$i"


done