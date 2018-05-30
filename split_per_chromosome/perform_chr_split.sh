#!/bin/bash

#Date 12/28/2017
#Take in .ped/.bed/.bim files and split into per chromosome 
#

usage() { #spit out the usage
cat <<UsageDisplay

perform_chr_split.sh -o <output dir> -i <input direct> 

Options:

-i input directory with .bim .bam. .fam file

-o output directory 

-x output prefix 

UsageDisplay

exit;

}

check_inputs()
{

    if [ "$out_pref" == "" ]
    then
	printf "Missing output prefix\n"
	usage;
    fi
    
    if [ "$output_dir" == "" ] || [ "$input_dir" == "" ] 
    then
	printf "Incorrect inputs\n"
	usage;
	
    fi
    
    if [ ! -d "$output_dir" ] || [ ! -d "$input_dir" ]
    then
	printf "Incorrect directory\n"
	usage;
    fi
    
    
    output_dir=$(cd $output_dir;pwd) #make full path ---
    input_dir=$(cd $input_dir;pwd)
    
}

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

out_pref=""
out_dir=""
input_dir=""

while getopts "x:i:o:h" args # iterate over arguments
do
    case $args in
	x)
	    out_pref=$OPTARG;; #prefix for output

	o)
	    output_dir=$OPTARG;;
	i)
	    input_dir=$OPTARG;;
	h)
	    usage;;
	*)
	    usage;;


    esac
done

check_inputs;

input_prefix=($(ls $input_dir/*.bim)) #get file name with bim
input_prefix=($(basename $input_prefix | rev | cut -c 5- | rev)) # get only the first name of *.bim file

for i in {1..22}
do

    out_dir=$output_dir/"CHR"$i
    mkdir $out_dir
    
    qsub -wd "$out_dir" -N "split_""$i" -o "$out_dir/"split_out_\$JOB_ID.out -e "$out_dir/"err_split\$JOB_ID.err -v bash /home/ss5505/scripts_cumc/split_per_chromosome/launch_job_split.sh -i $input_dir/$input_prefix -o $out_dir -c $i -x $out_pref

    
done
