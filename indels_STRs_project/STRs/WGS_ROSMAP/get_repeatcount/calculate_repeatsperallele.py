#!/usr/bin/env python

###Date 05/13/2019
#Sanjeev Sariya
###

##store repeat info from  bed file
#### python calculate_repeatsperallele.py
import re,sys,pprint

input_chr=10
input_bedfile="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/calculate_dosage.byrepeats/inputBED.repeats"

#get VCF info per CHR and read it. # zgrep -v  "#" merged_ROSmap_hipSTRCHR10.vcf.gz | cut -f 1-5 
input_STRalleles="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/calculate_dosage.byrepeats/inputtest_CHR10new"

##store STR name and their motifs in a dictionary

strNames_motif={} ### get details per STR name and their motif

with open(input_bedfile) as handle_bed:

    for line in handle_bed:
        
        line=line.rstrip()
        names_array=line.split()
        try:

            if len(names_array)!=7: ##hard coded length
                print "We've issues with ",line,len(names_array)
                next
            
        except Exception, err_arraySplit_length:
            print err_arraySplit_length
        
        ##if there are issues with length of array and splitting
        try:
            
            if int(names_array[3])!=len(names_array[6]):
                print "Issues with length of STR",line
                next
            ##if check motif length matches or not
            
        except Exception, err_motifLength:
            print err_motifLength,line
        ##if length has issue of motif and one stored in BED data
        
        strNames_motif[names_array[5]]=names_array[6]
    ##for loop ends
##with ends

print "we have stored STR info in dict"

###STR motif consolidation are over.

store_motifcounts={} ## key-value./ value is tab seprated for ref and allele1,allele2 
with open(input_STRalleles) as alleles_handle:
    for line in alleles_handle:
        
        line=line.rstrip()
        allelenames_array=line.split()
        
        if allelenames_array[4]==".":
            next
        ##if alt allele is not identified . Skip everything.

        else:
            ##if alt allele is identified
            
            str_name=allelenames_array[2]
            ref_allele=allelenames_array[3]
            alt_allele=allelenames_array[4]
            
            if str_name in strNames_motif:
                count_stringtostore=""
                
                motif=strNames_motif[str_name]
                ref_count=ref_allele.count(motif)
                count_stringtostore=str(ref_count) ##store as string
                
                ##now work on alternate alleles
                if ',' in alt_allele:
                    ##if multiple allele in alt 
                    array_multi_allele=alt_allele.split(',')

                    for val in range(0,len(array_multi_allele)):

                        if val==0:
                            temp_count=(array_multi_allele[val]).count(motif)
                            count_stringtostore=count_stringtostore+"\t"+str(temp_count)

                        else:
                            temp_count=(array_multi_allele[val]).count(motif)
                            count_stringtostore=count_stringtostore+","+str(temp_count)

                    #for loop ends for STRs with multi alleles
                    
                    store_motifcounts[str_name]=count_stringtostore
                else:
                    ##else if single alt allele
                    alt_count=alt_allele.count(motif)
                    count_stringtostore=count_stringtostore+"\t"+str(alt_count)
                    store_motifcounts[str_name]=count_stringtostore
                #check ends for , in alleles
            ##if check STR name exists in dict
                
            else:
                next
            ## if ST not find in dict check ends
        ##else if alternate allele is not .
            
    ##for loop ends
##with loop for alleles end

for key in store_motifcounts:
    print key,store_motifcounts[key]
    
