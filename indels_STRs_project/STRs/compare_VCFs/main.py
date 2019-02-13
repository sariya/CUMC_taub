#!/usr/bin/env python

#
#Date 02-05-2019
#Sanjeev Sariya


#Name: PyVCF
#Version: 0.6.8

#
#Maintain INFO string by recode-all-Info flag in vcftools. 
#Both these may or may not be multi-allelic. Depends on code and day you're doing. 
#This script assumes data are cleaned, filtered for calling rate, depth, MAF and other metrics necessay to
#perform analysis
#Takes two VFs - plink VCF from imputation and VCF from STR called.

#
#Python 2.7.13 (default, Jan 19 2017, 14:48:08)
#[GCC 6.3.0 20170118] on linux2
#

# fuzzywuzzy  Version: 0.17.0

import vcf,sys

def get_pos(temp_vcffile):
    
    """
    Read VCF file and send unique positions
    """
    
    pos_store=set() #store and return later
    vcf_reader_impute = vcf.Reader(filename=temp_vcffile)
    
    for record in vcf_reader_impute:
        pos_store.add(record.POS)

    #--for loop end

    return (pos_store)
    #---------------------
    #func ends
    #-----------------------    

if __name__=="__main__":

    imputed_vcf_file="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/IMPUTED_CHR22_batch2/plink_STRs_68_0.4.vcf"    
    called_file_vcf="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/bam_files_processing_STR/comparevcfs/updated_ids68.vcf"

    sitepos_imputation=get_pos(imputed_vcf_file) #get site position in a set
    print(len(sitepos_imputation))
    print "Sites for imputation completed"
    
    sitepos_vcalled=get_pos(called_file_vcf)
    print(len(sitepos_vcalled))
    print "Sites for seq completed"
    
    valuuuu=(sitepos_imputation).intersection(sitepos_vcalled)
    print len(valuuuu)

    print "Done reading"

##
