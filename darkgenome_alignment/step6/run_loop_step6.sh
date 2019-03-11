#!/bin/bash

#Date 03/07/2019
#Sanjeev Sariya

cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/sequencing_projects/darkgenome/

#
#Take expanded bed file and i from the other script.
#The i is multiplied by 2 for ploidy and later used in downstream processes in GATK for read mapping as stated in paper.
#

#$1 is expanded bed file and $2 is integer from 2 1 5 loop
REF="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/test_symlinks/human_g1k_v37_decoy.fasta"
RESULT_DIR="./TEST_DIR"
PREFIX="CR1.region"

genome="./genome_CHR_size.bedtools"
ploidy=$((2*$2))
out="${RESULT_DIR}/${PREFIX}.ploidy_${ploidy}.fa"

printf "$1 $2 $PREFIX $ploidy $out\n"

set -x 
awk "NF == $2" "$1" | /mnt/mfs/cluster/bin/bin/bedtools complement -i - -g $genome | /mnt/mfs/cluster/bin/bin/bedtools maskfasta -fi $REF -bed - -fo $out

bwa index -a bwtsw $out
set +x

printf "Completed bwa index \n"
set -x
samtools faidx $out
set +x
printf "Completed samtools index \n"

OUTPUT_dict=${out%.*}.dict #use this for picard

printf "$OUTPUT_dict\n"
set -x 
java  -Xmx20G -jar /mnt/mfs/cluster/bin/bin/picard.jar CreateSequenceDictionary  REFERENCE=$out  OUTPUT=$OUTPUT_dict
set +x
printf "Completed picard jar\n"

printf "Completed with $1 $2\n"

