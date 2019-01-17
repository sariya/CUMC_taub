#!/bin/python

#
#Date 01-17-2019
#Sanjeev Sariya
#
#Input for sets in GMMAT rv tests #we have weight as 1
#Output is a text file with prefix and path from user. This script is run per chromsome
#plink.annot has removed SNPs with . gene annotation
#python create_groups_snpsfiltered.py -i snps_annotation_VEP -d ./ -p batch1_CHR1

#
#Input file is
#21:31709731:G:A ENSG00000206107 KRTAP27-1
#21:31744127:A:T ENSG00000182816 KRTAP13-2
#21:31768494:C:A ENSG00000198390 KRTAP13-1
#21:35893717:G:T ENSG00000159200 RCAN1
#
#
#Make output as LOC646070       6       148277  C       T       1
#
#
" Python 2.7.13 (default, Jan 19 2017, 14:48:08) [GCC 6.3.0 20170118] on linux2 "
import argparse,os

if __name__=="__main__":
    print "Hello"
    parser=argparse.ArgumentParser("parse plink annoted files")
    
    parser.add_argument ('-i','--inp_annot',help='location/path of input annotated file',required=True) # store input plink annoated file
    parser.add_argument ('-d','--dir_out',help='location/path of output ',required=True) # store path for output dir
    parser.add_argument ('-p','--prefix_output',help='prefix for output',required=True) #
    args_dict = vars(parser.parse_args()) # make them dict..
    ####
    input_file=os.path.abspath(args_dict['inp_annot'])
    output_dir=os.path.abspath(args_dict['dir_out'])
    prefix_output=args_dict['prefix_output']
    output_file=output_dir+"/"+prefix_output+"groupfile.txt"
    with open(input_file) as handle:
        #
        for l in handle:
            line=l.strip() # strip new line
            space_array=line.split() # split line by space
            snp=space_array[0] #store SNP
            snp_split_array=snp.split(":") # split imputed SNP

            ##store values
            chromosome=snp_split_array[0]
            position=snp_split_array[1]
            ref_allele=snp_split_array[2]
            alt_allele=snp_split_array[3]
            gene_name=space_array[2]
            string_to_print=gene_name+"\t"+chromosome+"\t"+position+"\t"+ref_allele+"\t"+alt_allele+"\t"+"1"+"\n"
            
            with open(output_file, 'a') as the_file:
                the_file.write(string_to_print)
            #--with ends                                                                             
                                            
        #for loop for handle ends
    #with ends
                    
    
#--if main ends

