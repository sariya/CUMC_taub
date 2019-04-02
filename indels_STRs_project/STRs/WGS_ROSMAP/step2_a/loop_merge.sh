#!/bin/bash

cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/merged_STRs_hipSTR

#1,2,12,14
#{6,12,13,15}

##22, 10, 11 failed.
for i in {4,5,7,8,9,16,17,18}

do
qsub -cwd -N mergedCHR"$i" -b y bash run_merge.sh "$i"
done




