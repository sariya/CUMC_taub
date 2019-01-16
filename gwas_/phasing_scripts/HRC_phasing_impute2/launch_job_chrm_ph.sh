#!/bin/bash

#Date : 04/11/2018
#Place: Dr. GT
#CUMC 19th Floor

#1) split chromosome
#2) align
#3) phase
#4) Impute your data!!

#This script is being called from perform alignments
#-----------------------

usage() { #spit out the usage
cat <<UsageDisplay

launch_job_chrm_ph.sh

Options:

-b <bed file per chromsome>

-m <map file from the reference panel from the chromosome>

-f <hap file from reference panel for the chromosome>

-l <legend file from reference panel for the chromosome>

-s <sample file from reference panel for the chromosome>

-p <prefix CHR/chr and number>
UsageDisplay

exit;

}

#-------------------
if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

check_inputs(){
    if [ ! -f $bed_file".bed"  ]
    then
	printf "Inoccrect bed file\n"
	usage;
    fi
    
}
#------------------------

bed_file=""
map_files=""
hap_file="" 
legend_file="" 
sample_file=""
prefix="" #used for haps file later on. Phasing files have this information

#----
while getopts "p:f:b:m:l:s:h" args # iterate over arguments
do
    case $args in
	p)
	    prefix=$OPTARG;;
	
	b)
	    bed_file=$OPTARG;; #bed file
	m)
	    map_file=$OPTARG ;; #map file
	f)
	    hap_file=$OPTARG;; #hap file

	l)
	    legend_file=$OPTARG;; #legend
     
	s)
	    sample_file=$OPTARG;; #sample
	h)
	    usage;;
	*)
	    usage;;
    esac
done
#-------------
check_inputs;


set -x 
/home/ss5505/bin/shapeit \
    -B $bed_file \
    -M $map_file \
    --input-ref $hap_file $legend_file $sample_file \
    --duohmm \
    --output-graph $bed_file"phased"$prefix"-duohmm.graph" \
    --output-max $bed_file"phased"$prefix"-duohmm" \
    --thread 16 \
    --seed 123456789 \
    -W 5 
    
