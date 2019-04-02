#!/bin/bash

#Date 02/20/2019
#Sanjeev Sariya

bedfile_path="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/input_BEDfiles/chunked_bedfiles/"


##for i in {9,10,11,16,17,18,19,20}
##for i in {19..22}
##for i in {4,5,7}

##for i in {8,9,10,11,16,17,18}

for i in 11
do

chunk_files=($(ls $bedfile_path/hg19.hipstr_referenceCHR"$i"_chunk*))

for x in "${chunk_files[@]}"
do

chunk_name=($(basename $x) ) 

qsub -cwd -N CHRagain"$i" -l h_vmem=28G -v bash ./run_hipstr.sh "$i" "$x"
exit
done

done

