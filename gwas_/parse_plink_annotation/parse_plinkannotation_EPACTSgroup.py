#!/bin/python

#
#Date 01/15/2019
#Sanjeev Sariya
#Take output from plink annotation for gene
#Create group file for EPACTS
#Input for sets in GMMAT rv tests #we have weight as 1
#Output is a text file
#plink.annot has removed SNPs with . gene annotation
#
#python parse_plinkannotation_EPACTSgroup.py -i plink.annot -d ./ -p batch1_CHR1
#

import argparse,os
if __name__=="__main__":
    
    parser=argparse.ArgumentParser("parse plink annoted files")
    parser.add_argument ('-i','--inp_annot',help='location/path of input annotated file',required=True) # store input plink annoated file
    parser.add_argument ('-d','--dir_out',help='location/path of output ',required=True) # store path for output dir
    parser.add_argument ('-p','--prefix_output',help='prefix for output',required=True) #
    args_dict = vars(parser.parse_args()) # make them dict..
    
    input_file=os.path.abspath(args_dict['inp_annot'])
    output_dir=os.path.abspath(args_dict['dir_out'])
    prefix_output=args_dict['prefix_output']

    gene_snp={} ##use this dict later on
    
    with open(input_file) as handle:
        #{with starts
        for l in handle:
            #{ for loop for handle starts
            line=l.strip() # strip new line
            if ":" in line:
                
                space_array=line.split() # split line by space
                snp=space_array[1] #store SNP
                snp_split_array=snp.split(":") # split imputed SNP
                #--store different variables
                
                chromosome=snp_split_array[0]
                position=snp_split_array[1]
                ref_allele=snp_split_array[2]
                alt_allele=snp_split_array[3]
                
                if "|" in line:
                    #{ if line has |
                    
                    gene_annotated=space_array[3]
                    split_gene_array=gene_annotated.split("|")
                    
                    for i in split_gene_array:
                        
                        index=i.find('(') # find index of (
                        gene_name=(i[0:index]).strip() #strip last space
                        #--make a string and print to out file
                        string_to_print=gene_name+"\t"+chromosome+":"+position+"_"+ref_allele+"/"+alt_allele
                        
                        ##
                        ##Check gene exists
                        ##
                        print string_to_print
                        ##with open(output_file, 'a') as the_file:
                        ##  the_file.write(string_to_print)
                    else:
                        #{ if line doesn't have |
                        
                        gene_annotated=space_array[3]
                        index=gene_annotated.find('(')
                        gene_name=(gene_annotated[0:index]).strip() #strip last space
                        
                        string_to_print=gene_name+"\t"+chromosome+"\t"+position+"\t"+ref_allele+"\t"+alt_allele+"\t"+"1"+"\n"
                        #with open(output_file, 'a') as the_file:
                        #   the_file.write(string_to_print)
                #----if check ends for |
                        
            #--if : ends
        
                        
