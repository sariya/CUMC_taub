#!/bin/bash

cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/merged_STRs_hipSTR/
printf "CHR$1\n"


###/mnt/mfs/cluster/bin/bin/vcf-concat /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/input_bam_symlinks/rosmapCHR"$1"_chunk*.vcf.gz | gzip -c > merged_ROSmap_hipSTRCHR"$1".vcf.gz

/mnt/mfs/cluster/bin/bin/vcf-concat /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/failedCHR_STR_calling/rosmapCHR"$1"_chunk*.vcf.gz  | gzip -c > merged_ROSmap_hipSTRCHR"$1".vcf.gz

printf "completed\n"