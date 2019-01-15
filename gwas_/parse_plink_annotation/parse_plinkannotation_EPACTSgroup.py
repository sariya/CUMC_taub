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
if __name__=="__main__":
    parser=argparse.ArgumentParser("parse plink annoted files")
    parser.add_argument ('-i','--inp_annot',help='location/path of input annotated file',required=True) # store input plink annoated file
    parser.add_argument ('-d','--dir_out',help='location/path of output ',required=True) # store path for output dir
    parser.add_argument ('-p','--prefix_output',help='prefix for output',required=True) #
    args_dict = vars(parser.parse_args()) # make them dict..

    
