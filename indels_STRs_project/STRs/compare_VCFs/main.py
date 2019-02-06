#!/usr/bin/env python


#
#Date 02-05-2019
#Sanjeev Sariya
#

#
#This script somehow doesn't give INFO string in VCF. Weird. :(
#Both these VCFs are normalized for multi-alleles
#This script assumes data are cleaned, filtered for calling rate, depth, MAF and other metrics necessay to
#perform analysis
#Takes two VFs - plink VCF from imputation and VCF from STR called.
#

import vcf

if __name__=="__main__":

    imputed_vcf_file="/home/ss5505/scripts_cumc/indels_STRs_project/STRs/compare_VCFs/plink_updated_ids_STRs_68_0.8.vcf"    
    called_file_vcf="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/bam_files_processing_STR/compare_VCFs/normalized_CHR22_STR_fixedVersion.vcf.gz"
    
    vcf_reader = vcf.Reader(filename=file_vcf)

    for record in vcf_reader:
        print record, record.start,record.end
        print record.num_hom_ref, record.num_het, record.num_hom_alt       
        break

#    file_vcf="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/bam_files_processing_STR/str_calls_38Fams_toolStutter.vcf.gz"
#    vcf_reader = vcf.Reader(filename=file_vcf)

    print "Done reading"

##
