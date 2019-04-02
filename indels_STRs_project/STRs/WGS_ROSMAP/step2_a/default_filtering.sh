#!/bin/bash


for i in {1..9}

	 ## 12..21
do

    printf "CHR$i\n"


/usr/bin/python \
	/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/tools_STRs/hipSTR/HipSTR/scripts/filter_vcf.py \
	--vcf \
		    merged_ROSmap_hipSTRCHR"$i".vcf.gz \
                              --min-call-qual         0.9 \
                              --max-call-flank-indel  0.15 \
                              --max-call-stutter      0.15  \
		      --min-call-allele-bias  -2 \
			      --min-call-strand-bias  -2 > filtered_CHR"$i".vcf
done
