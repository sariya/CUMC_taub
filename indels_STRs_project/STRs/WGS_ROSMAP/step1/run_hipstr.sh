#!/bin/bash

#Date 02/20/2019
#Sanjeev Sariya

ulimit -S -n 4096

printf "CHR received $1\n"
printf "bed file received $2\n"

basename_bedfile=($(basename "$2"))
temp_namechunk=($(echo $basename_bedfile | sed 's/\.bed//g' | sed 's/hg19.hipstr_reference//g')) #use this in output vcf str called output
printf "chunk ends with $temp_namechunk\n"

inputfile="input_bampath.txt"

outputfile=rosmap"$temp_namechunk".vcf.gz 
printf "$outputfile\n"

#/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/tools_STRs/hipSTR/HipSTR/HipSTR --chrom "$1" --output-filters --bams "$temp_string"  \
#--fasta  /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/test_symlinks/human_g1k_v37_decoy.fasta #--min-reads 100 --regions "$2" \#--str-vcf $outputfile

set -x 

/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/tools_STRs/hipSTR/HipSTR/HipSTR \
--chrom "$1" --output-filters  \
--bam-files "kill_new.txt" \
--fasta /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/test_symlinks/human_g1k_v37_decoy.fasta \
--min-reads 100 \
--regions "$2" \
--str-vcf $outputfile

printf "Complete\n"


