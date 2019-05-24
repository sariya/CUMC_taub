#!/bin/bash

#
#Date 03/09/2019
#Sanjeev Sariya ROS MAP STR calls get basic stats
#

cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/merged_STRs_hipSTR

#{1,2,12,14}
# {6,12,13,15}
##{19,20,21,22}
for i in  {11,22}
do
    printf "CHR$i\n"
    qsub -cwd -N statsCHR"$i" -b y bash ./perform_basicdepth.sh "$i"
done

