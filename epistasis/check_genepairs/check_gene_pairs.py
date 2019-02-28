#!/usr/bin/env python

"""
Date 02/28/2019
Take in metal files an many provided by used and check for pairs such as:

PTK2B-CASS4
CASS4-PTK2B
Take in file with per line per file - it should be its PATH
We will take gene list from use 
Provide full path for file per line


"""

import os, argparse,sys
from itertools import combinations
import pprint

def getgene_combinations(tempgeneFile):
    """
    return dict for gene pairs. Key is number 1..N 
    where N is total number of gene pairs
    """
    gene_list=[] #store gene names in this and use for combinations
    pairs_genes={}
    with open(tempgeneFile) as handle:
        for line in handle:
            line=line.strip()
            gene_list.append(line)
        #for loop ends
    #with loop ends

    ##use gene list to create combinations
    
    #https://stackoverflow.com/a/17434095/2740831
    count=0 #store this in key for gene pairs
    
    for group in combinations(gene_list, 2):
        ##print('-'.join(group))         ##print(group)         print group[0],group[1]
        """
        group is tuple stored. We'll keep group[0] as Gene1 and group[1] as gene2
        gene1-gene2 combination. Assuming the tuple stays as is.

        """
        count+=1
        pairs_genes[count]=group #pprint.pprint(dictionary)
    #for loop ends

    return pairs_genes

##function ends----------------------------------------------------->>

def maxindex_(temp_string):

    return max([i for i, a in enumerate("HLA-DRB1-CLU") if a == "-"])

    ##Function ends
    ##
    
def minindex_():
    print ""
    ##
    ##
def read_file(temp_filelist):
    
    """
    iterate through file list.
    Open each file and read. Then look for pairs
    """
    for i in temp_filelist:
        
        print i
        with open(i) as handle:
            for line in handle:
                line=line.strip()
                infile_pair=(line.split())[0]
                print infile_pair
                """
                For the gene pair now look for dash, and see where to split. 
                Count # of dashes. to make a cut
                """
                try:
                    
                    count_dash=infile_pair.count('-')
                    gene1="" #split by dash and store gene1
                    gene2="" #split by dash and gene2
                    
                    if count_dash==1:
                        print "Two paired non HLA gene ",infile_pair
                        array_split_gene=infile_pair.split('-')
                        gene1=array_split_gene[0]
                        gene2=array_split_gene[1]
                        
                    elif count_dash==2:
                        print "HLA gene with some unique gene ",infile_pair

                        if infile_pair.startswith("HLA"):
                            print "we have HLA in beginning",infile_pair.upper()
                        else:
                            if "HLA" in infile_pair:
                                print "we have something different",line

                            else:
                                print "Something's wrong. No beginning with HLA",line
                    elif count_dash==3:
                        print " HLA gene with HLA gene ",infile_pair
                    else:
                        print "We don't know what's up here!",line
                        sys.exit()
                    ##if elif ends
                    
                ##try ends
                except Exception as error_dash:
                    print error_dash
                ##except ends

                
            ##for loop ends for iterating over file
        ##with loop ends
        print "Completed file ",
    #for loop ends for file iterations

    ##def ends for read file 
if __name__=="__main__":

    parser=argparse.ArgumentParser("get gene pairs")
    parser.add_argument('-g','--gene',required=True,type=str)
    parser.add_argument('-m','--metal',required=True,type=str)
    args_dict = vars(parser.parse_args()) # make them dict..
    metal_listfile= args_dict['metal']
    gene_file=args_dict['gene']
        
    if not os.path.isfile(gene_file):
        print "Missing files"
        sys.exit()
    #--if check ends
    gene_pairs=getgene_combinations(gene_file)
    pprint.pprint(gene_pairs)
    file_list=[] #store name of files in this . lists are ordered
    with open(metal_listfile) as handle:
        for line in handle:
            temp_file= line.strip()
            
            if not os.path.isfile(temp_file):
                print "Missing files"
                sys.exit()
            else:
                file_list.append(temp_file)
            #if for file check
            
        #for loop ends
    #with loop ends
    print len(file_list)," The number of files are "
    read_file(file_list)
