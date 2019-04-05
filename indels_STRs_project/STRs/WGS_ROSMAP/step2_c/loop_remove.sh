#!/bin/bash

cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/cleaning_mergedSTRs/

for i in 10
do
printf "CHR$i\n"

##bash ./perform_removeindividuals.sh "$i"

qsub -cwd -N removechr"$i" -v bash ./perform_removeindividuals.sh "$i"

done
