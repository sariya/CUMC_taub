#!/usr/bin/env python
##Date 04/02/2019
##Sanjeev Sariya
##Python 2.7.13 (default, Jan 19 2017, 14:48:08)
##[GCC 6.3.0 20170118] on linux2

import vcf,sys,re
import numpy as np
vcf_file="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/cleaning_mergedSTRs/filter.str.no_monomorphic.poor_depth.highmissing.postcleaning.removeindi_CHR21.recode.vcf"

vcf_reader = vcf.Reader(open(vcf_file, 'r'))
record = next(vcf_reader)
##/home/ss5505/scripts_cumc/indels_STRs_project/STRs/WGS_ROSMAP/step3

row_counts=0 
for information in vcf_reader:
    row_counts+=1
##

print "Found total records as ",row_counts

sample_dict={} ##store samples in dictionary

itr=1 ##use this for interation and storing.
for x in vcf_reader.samples:
    sample_dict[itr]=x
    itr+=1
    
#for loop ends

itr=1


#record = next(vcf_reader)
##https://buildmedia.readthedocs.org/media/pdf/pyvcf/latest/pyvcf.pdf

store_genotype=np.zeros((row_counts,len(sample_dict))) ##create empty matrix
store_variantID=np.chararray((row_counts,3)) ##create empty matrix chararray

itr=1 ##iterate though samples 
itr_row=0
vcf_reader = vcf.Reader(open(vcf_file, 'r'))

for temp in vcf_reader:
    print temp
    
    store_variantID[itr_row,0]=temp.ID+":"+str(temp.CHROM)+":"+str(temp.POS)
    print temp.ID+":"+str(temp.CHROM)+":"+str(temp.POS)
    store_variantID[itr_row,1]="YX"
    store_variantID[itr_row,2]="XY"
    print len(store_variantID[itr_row,0])

    sys.exit()
    itr_sample=0
    for sample in temp.samples:

        if sample['GT'] =="./." or  sample['GT']=="0|0" or sample['GT'] ==".":
        ##we don't do anything and keep genotype as zeros
            pass
    
        else:
            split_genotype=(sample['GT']).split("|") #split genotype check ref genotype first 

            ##make first element as interger. Convert String into integer
            if int(split_genotype[0]) == 0:
                store_genotype[itr_row,itr_sample]=1 #if first allele is 0 then dominant one genotype
                
            if  int(split_genotype[0]) !=0:
                store_genotype[itr_row,itr_sample]=2


        itr_sample+=1
            #if ends
        ##else ends

    ### for loop ends for samples
    print "Total samples are ",itr_sample
    print store_genotype.shape

    temp_string="" ##
    
    print store_variantID[itr_row,]
    for o in store_genotype[itr_row,]:
        temp_string=temp_string+","+str(o)
    print temp_string
    itr_row+=1
    sys.exit()

###

for sample in record.samples:
    
    if sample['GT'] =="./." or  sample['GT']=="0|0" \
       or sample['GT'] ==".":
        ##we don't do anything and keep genotype as zeros
        pass
    
    else:
        split_genotype=(sample['GT']).split("|") #split genotype check ref genotype first 

        if split_genotype[0] == 0:
            store_genotype[itr_row,itr]=1 #if first allele is 0 then dominant one genotype
            print "We have something new",sample['GT'],sample,sample_dict[itr]
            print itr_row,itr
            
        if  split_genotype[0] !=0:
            store_genotype[itr_row,itr]=2
            print "We have something new",sample['GT'],sample,sample_dict[itr]
            print itr_row,itr
        #if ends
    ##else ends
        
    itr_row+=1
    itr+=1
#for loop ends

print "Evrything is complete!!!"

    
