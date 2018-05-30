#!/bin/bash

#Date Nov 27 2017
#Sanjeev Sariya
#Dr. Giuseppe Tosto


#
#Output will be made in the Chromosome folders

#
check_inputs(){
if [ "$input_dir" == "" ] || [ "$chrom_pre" == "" ] || [ "$prefix_input_files" == "" ] || [ "$ref_path" == "" ]
then
    printf "Incorrect inputs\n"
    usage;
    
fi

if [ ! -d "$input_dir" ] || [ ! -d "$ref_path" ]
then 
    printf "Incorrect directory"
    usage;
fi


input_dir=$(cd $input_dir;pwd)
ref_path=$(cd $ref_path;pwd)
}

usage() { #spit out the usage
cat <<UsageDisplay

perform_alignment.sh 

Options:

-i <path of input directory where all chromosome folders are present>

-c <chr|CHR> Folder names start with either CHR or chr

-p <prefix of .bim .bed file> 

-r <path for reference files>

bash ~/scripts_cumc/gwas_/alignment/perform_alignment.sh -i /home/ss5505/tool_toy_files/shapeit_impute2/test_data -r ~/reference_panels/1000GP_Phase3 -c CHR -p HGWAS6_final_

UsageDisplay

exit;

}

if [ $# -eq 0 ] # exit if no arguments!
then
    usage;
fi

prefix_input_files="" #store name of .bim/.bed/.fam file
chrom_pre=""  #chr or CHR or CHROM
input_dir="" #path where chr folders are present
ref_path="" #direcotry where all reference files are present

while getopts "r:i:c:p:h" args # iterate over arguments
do
    case $args in
	r)
	    ref_path=$OPTARG;;
	i)
	    input_dir=$OPTARG;;
	p)
	    prefix_input_files=$OPTARG;;
	c)
	    chrom_pre=$OPTARG;;
	h)
	    usage;;
	*)
	    usage;;
    esac
done

check_inputs;

for chr_num in {1..22}
do

    #--loop through chromosomes
    
    #--reference files

    #--set bed file per chromsome and get hap, legend, map, sample files from reference panels
    
    b_file=$input_dir/$chrom_pre$chr_num/$prefix_input_files$chrom_pre$chr_num #store bed/bim/fam file

    m_file=$ref_path/genetic_map_chr${chr_num}_combined_b37.txt  #store map file
    hap_file=$ref_path/1000GP_Phase3_chr${chr_num}.hap.gz #store hap file
    legend_file=$ref_path/1000GP_Phase3_chr${chr_num}.legend.gz #store legend file
    sample_file=$ref_path/"1000GP_Phase3.sample" #store sample file

    #----

    output_log=$input_dir/$chrom_pre$chr_num/$chrom_pre$chr_num".alignments" #output name for alignments logs

    #--submit job 
set -x     
    qsub -wd "$input_dir/$chrom_pre$chr_num/" -o "$input_dir/$chrom_pre$chr_num/"algn_out_\$JOB_ID.out -e "$input_dir/$chrom_pre$chr_num/"algn_err_\$JOB_ID.err -N "algn_chr""$chr_num" -v bash /mnt/mfs/hgrcgrid/homes/ss5505/scripts_cumc/gwas_/alignment/1KGP/launch_job_chrom.sh -b $b_file -m $m_file -f $hap_file -l $legend_file -s $sample_file
    
done
