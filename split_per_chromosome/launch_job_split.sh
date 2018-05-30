#!/bin/bash

source ~/.bashrc

usage() { #spit out the usage
cat <<UsageDisplay

launch_job_split.sh -o <output dir> -i <input direct> -c <chromsome>

Options:

-i <input_prefix>

UsageDisplay

exit;

}

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi


while getopts "x:c:i:o:h" args # iterate over arguments
do
    case $args in
	c)
	    chr_name=$OPTARG;;
	x)
	    prefix_out=$OPTARG;; #use it for output prefix
	o)
	    out_dir=$OPTARG;; #directory for output
	i)
	    infile_base=$OPTARG;;

    esac
done

cd $out_dir

plink --nonfounders --allow-no-sex --bfile $infile_base --chr $chr_name --make-bed --out "$prefix_out""_CHR""$chr_name"  
