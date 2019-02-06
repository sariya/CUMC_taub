#!/usr/bin/env python


#
#Date 02-05-2019
#Sanjeev Sariya
#

import vcf

if __name__=="__main__":
    print "Hello\n"

    file_vcf="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/bam_files_processing_STR/compare_VCFs/normalized_CHR22_STR_fixedVersion.vcf.gz"
    vcf_reader = vcf.Reader(filename=file_vcf)


    for record in vcf_reader:
        print record, record.INFO, record.FORMAT,vars(record)
        print record.num_hom_ref, record.num_het, record.num_hom_alt
        
        break

    file_vcf="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/bam_files_processing_STR/str_calls_38Fams_toolStutter.vcf.gz"
    vcf_reader = vcf.Reader(filename=file_vcf)
    for record in vcf_reader:
        print record, record.INFO,record.INFO['DFLANKINDEL']
        print record.num_hom_ref, record.num_het, record.num_hom_alt
        
        break
    print "Done reading"

##
