#!/usr/bin/env python

#
#Date 02-05-2019
#Sanjeev Sariya


#Name: PyVCF
#Version: 0.6.8

#Imptuation data are reference and for which comparisons are made
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

def find_only_calledsites(set_impute_sites, set_strcalled_sites):

    print len(set_impute_sites)
    print len(set_strcalled_sites)

    """
    We'll take sites from VCF and take sites - imputed and STR called. 
    We'll find sites that are only present in STR called data
    
    """
    
    #work on element in impute but not in VCF called https://stackoverflow.com/a/50559206/2740831
    sites_calledonly=(set_strcalled_sites^set_impute_sites)&set_strcalled_sites
    print(len(sites_calledonly))

    calledSTR_calledonly_list = list(sites_calledonly) #convert set into list to access elements # https://stackoverflow.com/a/32512277/2740831

    print calledSTR_calledonly_list[1:10]
    print "Leaving only called STRs function"
    
    #
    #function ends
    #<<->>
    
def find_only_imputationsites(set_impute_sites, set_strcalled_sites):

    print len(set_impute_sites)
    print len(set_strcalled_sites)

    """
    We'll take sites from VCF and take sites - imputed and STR called. 
    We'll find sites that are only present in Imputed data
    
    """
    
    #work on element in impute but not in VCF called https://stackoverflow.com/a/50559206/2740831
    sites_imputedonly=(set_impute_sites^set_strcalled_sites)&set_impute_sites
    print(len(sites_imputedonly))

    imputed_sitesonly_list = list(sites_imputedonly) #convert set into list to access elements # https://stackoverflow.com/a/32512277/2740831
    #print sites_imputedonly[1:10]

    print imputed_sitesonly_list[1:10]
    print "Leaving only impute SNPs function"

    #
    #function ends
    #<<->>
    
if __name__=="__main__":

    imputed_vcf_file="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/IMPUTED_CHR22_batch2/plink_STRs_68_0.4.vcf"    
    called_file_vcf="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/bam_files_processing_STR/comparevcfs/updated_ids68.vcf"

    sitepos_imputation=get_pos(imputed_vcf_file) #get site position in a set
    print(len(sitepos_imputation))
    print "Sites for imputation completed"
    
    sitepos_vcalled=get_pos(called_file_vcf)
    print(len(sitepos_vcalled))
    print "Sites for seq completed"
    
    overlapped_sites=(sitepos_imputation).intersection(sitepos_vcalled)
    my_overlapped_list = list(overlapped_sites) #convert set into list to access elements # https://stackoverflow.com/a/32512277/2740831
    print len(my_overlapped_list)
    print my_overlapped_list[1:5]

    #find sites unique to STR called and imputation

    find_only_imputationsites(sitepos_imputation,sitepos_vcalled) # sent to find sites only imputation
    find_only_calledsites(sitepos_imputation,sitepos_vcalled)

    #
    print "Done reading"

    

##
