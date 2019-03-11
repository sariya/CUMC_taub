#!/bin/bash

#
#Date 03/06/2019
#Sanjeev Sariya

cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/sequencing_projects/darkgenome/

camo_bed="./CR1_gene.bed" #input bed file.

genome="./genome_CHR_size.bedtools" #this is used as bed genomic boundaries . 
expanded="CR1_gene_expanded.bed" #use it for exapnded boundaries 

#expand input bed file 
/mnt/mfs/cluster/bin/bin/bedtools slop -b 50 -i $camo_bed -g $genome > $expanded

#
#the 5th column in shared bed file tells number of regions boundary spans. 2, 3,4 based on this number the ploidy is decided. 
#say if region contains 3 exon/introns. Ploidy is 6. That is 1/6 of the reads must map to the region 
#

for i in $(seq 2 1 5)
do

qsub -cwd -N ploidy"$i" -b y bash run_loop_step6.sh $expanded "$i"

done

exit
