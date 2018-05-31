#!//usr/bin/python
__date__="May 30 2018"
__location__="CUMC TaubNY, NY"
__author__="Sanjeev Sariya | Dr. Tosto"
__version__="0.1"
"""

Take output from AKT. Make sure to print properly. Parse _
Produce an output that's tab delimited (transposed)
python parse_duplicate.py -a allrelatives  > again_dup_check

"""

import os, argparse,re,sys

#---------------------------------------------------------------->>>
if __name__=="__main__":

	parser=argparse.ArgumentParser("description")
        parser.add_argument('-a','--akt',required=True,type=str)
	args_dict = vars(parser.parse_args()) # make them dict..
        akt_file= args_dict['akt']
        
        dupl_dict={} #store dups in this dict
        #--print this in the ends
        with open(akt_file) as handle:
                for line in handle:

                        line=line.rstrip()
                        dup_names=re.split('\s+',line)

                        if "Dup" in line:
                                #---if dup in line --then proceed ahead
                                
                                count_under=line.count('_')
                                dup_names=re.split('\s+',line)
                                iid_underscore=re.split('_',dup_names[1])
                                dup_string="" #store for dict later wards

                                if len(iid_underscore)==2:
                                        
                                        #--if only one underscore
                                        if iid_underscore[0] == "0":
                                                dup_string=iid_underscore[1]
                                        else:
                                                dup_string=iid_underscore[0]+"_"+iid_underscore[1]                                                
                                #--check for one underscore ends
                                
                                if len(iid_underscore)>2:
                                        dup_string=iid_underscore[1]+"_"+iid_underscore[2]

                                #----check fore more than 2 _ ends

                                #--check key existance in dict. 
                                if dup_names[0] in dupl_dict:   
                                        dupl_dict[dup_names[0]]=dupl_dict[dup_names[0]]+"\t"+dup_string
                                else:   
                                        dupl_dict[dup_names[0]]=dup_string
                                #--store in dict                                        
                        #--if else ends for Key in dictionary
                #--for loop ends ---<>>>>>>>>>>>>>>>>>>>>
        #--with ends

        #--print dict
        
        for key in dupl_dict:
                print key+"\t"+dupl_dict[key]
        #--for print ends
