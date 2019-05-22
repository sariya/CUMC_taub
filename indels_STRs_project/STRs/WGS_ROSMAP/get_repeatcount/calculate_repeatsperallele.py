#!/usr/bin/env python

###Date 05/13/2019
#Sanjeev Sariya
###

##store repeat info from  bed file
#### python calculate_repeatsperallele.py
import re,sys, os, argparse
from pprint import pprint
from collections import OrderedDict

parser=argparse.ArgumentParser("Get motif repeat counts in alleles")

parser.add_argument ('-b','--bed',help='location/path of bed file',required=True) # input bed file 
parser.add_argument ('-s','--str',help='location/path of str file',required=True) # str file
parser.add_argument ('-d','--dir_out',help='location/path of output ',required=True) # store path for output dir
parser.add_argument ('-p','--prefix_output',help='prefix for output',required=True) # provided by loop script

args_dict = vars(parser.parse_args()) # make them dict..
                    

#output_file="testcomplete"
#input_bedfile="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/input_BEDfiles/unchunked_bedfiles/hg19.hipstr_reference_CHR18.bed"

#get VCF info per CHR and read it. # zgrep -v  "#" merged_ROSmap_hipSTRCHR10.vcf.gz | cut -f 1-5 
#input_STRalleles="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/calculate_dosage.byrepeats/input_getcount_CHR18_STR.vcfdata"

input_bedfile=os.path.abspath(args_dict['bed'])
output_file=os.path.abspath(args_dict['dir_out'])+"/"+args_dict['prefix_output']
input_STRalleles=os.path.abspath(args_dict['str'])
                        
print input_bedfile,output_file,input_STRalleles

#output_file="test" #input_chr=10
##input_bedfile="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/calculate_dosage.byrepeats/temp"

##store STR name and their motifs in a dictionary
strNames_motif=OrderedDict() ### get details per STR name and their motif

with open(input_bedfile) as handle_bed:

    for line in handle_bed:
        
        line=line.rstrip()
        names_array=line.split()
        try:

            if len(names_array)!=7: ##hard coded length
                print "We've issues with ",line,len(names_array)
                continue
            
        except Exception, err_arraySplit_length:
            print err_arraySplit_length
        
        ##if there are issues with length of array and splitting
        try:

            if names_array[6]=="N/A":
                print "There's nasty stuff here", line
                continue
            ##if check motif length matches or not
            
            elif names_array[3]=='1':
                print "motif of length one ",line
                continue
            
            elif int(names_array[3])!=len(names_array[6]):
                strNames_motif[names_array[5]]=names_array[6]
            
            else:
                strNames_motif[names_array[5]]=names_array[6]

        except Exception, err_motifLength:
            print err_motifLength,line
        ##if length has issue of motif and one stored in BED data
        
    ##for loop ends
##with ends

print "we have stored STR info in dict"
print len(strNames_motif)
##pprint(strNames_motif )

###STR motif consolidation are over.

store_motifcounts=OrderedDict() ## key-value./ value is tab seprated for ref and allele1,allele2

################
with open(input_STRalleles) as alleles_handle:
    for line in alleles_handle:
        
        line=line.rstrip()
        allelenames_array=line.split()
        
        if allelenames_array[4]==".":
            continue
        ##if alt allele is not identified . Skip everything.

        else:
            ##if alt allele is identified
            
            str_name=allelenames_array[2]
            ref_allele=allelenames_array[3]
            alt_allele=allelenames_array[4]
            
            if str_name in strNames_motif:
                count_stringtostore=""


                if "/" in strNames_motif[str_name]:
                    
                    ##work with multiple motifs
                    
                    multi_motif=strNames_motif[str_name]
                    
                    ##split by /
                    array_multiMotif=multi_motif.split('/')
                    ##get counts in ref allele first
                    counts_ref=0
                    
                    for x_motif in array_multiMotif:
                        counts_ref+=ref_allele.count(x_motif)
                        
                    ##for loop ends for getting multi motifs in  ref allele
                    
                    count_stringtostore=str(counts_ref)

                    if ',' in alt_allele:
                        
                        ##if multi alt alleles
                        array_multi_allele=alt_allele.split(',')
                        
                        for val in range(0,len(array_multi_allele)):
                            ##iterate over multi allele
                            
                            counts_alt=0
                            
                            if val==0:
                                
                                for x_motif in array_multiMotif:
                                    counts_alt+=(array_multi_allele[val]).count(x_motif)
                                #for loop ends for motif iteration
                                
                                count_stringtostore=count_stringtostore+"\t"+str(counts_alt)
                            else:
                                
                                for x_motif in array_multiMotif:
                                    counts_alt+=(array_multi_allele[val]).count(x_motif)
                                    
                                ##for loop ends for motif iteration
                                count_stringtostore=count_stringtostore+","+str(temp_count)
                            ##if val is not 0
                            
                        #for loop ends for STRs with multi alleles
                    ###if check ends for multi alt alleles
                    
                    else:

                        ##if only one alt allele
                        counts_alt=0
                        
                        for x_motif in array_multiMotif:
                            counts_alt+=alt_allele.count(x_motif)
                            
                        ##for loop ends for counting multi motifs in alt allele
                        
                        count_stringtostore=count_stringtostore+"\t"+str(counts_alt)
                        store_motifcounts[str_name]=count_stringtostore
                                                                        
                    #####we work with multiple motifs
                #########################################################################################################
                else:
                    
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
                print "STR not found ", line
                continue
            ## if ST not find in dict check ends
        ##else if alternate allele is not .
            
    ##for loop ends
##with loop for alleles end

#write to output file.

print len(store_motifcounts)

with open(output_file, 'a') as handle_write:    
    for key in store_motifcounts:
        handle_write.write(key+"\t"+store_motifcounts[key]+'\n')
    #for loop ends
##with loop ends

print "completed looping, writing and calculating counts"
