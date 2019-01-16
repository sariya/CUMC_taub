#!/bin/bash

###
###Date 03 07 2018
###
input_dir=""
merge_file=""
out_dir=""
chr=""

while getopts "c:o:m:i:h" args # iterate over arguments
do
    case $args in
	c)
	    chr=$OPTARG;; #store chromosome
	
	m)
	    merge_file=$OPTARG;; #output merged file name

	i)
	    input_dir=$OPTARG;; #input dosage dir
	o)
	    out_dir=$OPTARG;; #output dir

	h)
	    usage;;
	*)
	    usage;;
    esac
done

printf "$input_dir\n$merge_file\n$out_dir\n"
temp_file=$out_dir/"temp_merge_chr$chr.text"
printf "$temp_file\n"

impute_info_file=( $(ls  -v $input_dir/*.impute2_info ) )

###
### Loop through the files
###
set -x

for f in "${impute_info_file[@]}"
do

    cat $f >> $temp_file
done

awk '!a[$0]++' $temp_file >> $merge_file
#rm $temp_file

printf "All done for $chr\n"
