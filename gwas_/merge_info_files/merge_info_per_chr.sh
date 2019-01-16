#!/bin/bash


# __Date__  Dec 06 2017
# __location__ CUMC 19th Floor
# __PI__ Dr. Giuseppe Tosto

#
#merge info Impute files  per chromosome of different chunks
#

usage(){
    
cat <<UsageDisplay

merge_info_per_chr.sh -i <input directory> -c <chr|chrom|CHR|CHROM> -o <output dir>

Options:

-i input directory where all CHR folder are present

-c prefix for chromosome individual files chr, CHROM, etc. based on plink data formatted or created folder. 

-o output dir - make this before submitting script

UsageDisplay

exit;

}

check_inputs()
{

if [ "$chr_prefix" == "" ]
then 
    printf "Chr prfix is missing\n"
    usage;
fi

if [ "$input_dir" == "" ] || [  "$chr_prefix" == "" ] || [ "$out_dir" == "" ] 
then
    printf "Incorrect inputs\n"
    usage;
    
fi

if  [ ! -d "$input_dir" ] || [ ! -d "$out_dir" ]
then
    printf "Incorrect directory\n"
    usage;
fi

input_dir=$(cd $input_dir;pwd)
out_dir=$(cd $out_dir;pwd)

}
#--input check ends ----

output_dir="" #store output dir
chr_prefix="" #generally files are outputted based on chr! we need that prefix.
out_dir=""

while getopts "o:c:i:h" args # iterate over arguments
do
    case $args in
	o)
	    out_dir=$OPTARG;;  
	c)
	    chr_prefix=$OPTARG;;
	
	i)
	    input_dir=$OPTARG;;

	h)
	    usage;;
	*)
	    usage;;
    esac
done

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

check_inputs;

for i in {1..22}
do
    chr_dir=$chr_prefix$i # get chr direcotrpy

    if [ -d $input_dir/$chr_dir ] 
    then

	dosage_dir=$input_dir/$chr_dir/"DOSAGE_FILES"

	if [ -d $dosage_dir ] 
	then
	    
	    impute_info_file=( $(ls  -v $dosage_dir/*.impute2_info ) )
	    
	    if [[ "${#impute_file[@]}" -ge 0 ]]
	    then

		merge_file_name=$out_dir/"merged_info_$chr_dir".info
#		bash /home/ss5505/scripts_cumc/gwas_/merge_info_files/launch_job.sh -m $merge_file_name -i $dosage_dir -o $out_dir -c $i 

		qsub -N "chr$i" -l h_vmem=15G -wd $out_dir -v bash /home/ss5505/scripts_cumc/gwas_/merge_info_files/launch_job.sh -i $dosage_dir -m $merge_file_name -o $out_dir -c $i

		
#		for f in "${impute_info_file[@]}"#		do#		    cat $f >> $merge_file_name#		done

	    else
		printf "$dosage_dir has issues with impute files \n"
	    fi

	else
	    printf "$chr_dir has issues with Dosage directory\n"

	fi
    else
    printf "$chr_dir isn't found in input directory\n"
    fi
    printf "Done through CHR$i\n"
done

printf "\n<--Complete-->\n"
