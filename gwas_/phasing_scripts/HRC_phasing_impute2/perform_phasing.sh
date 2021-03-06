#!/bin/bash

#Date April 11  2018
#Sanjeev Sariya
#Dr. Giuseppe Tosto


#
#Output will be made in the Chromosome folders  we need to have phasigin using HRC panel
##CHR 1 is a different animal here. with No IBD data
##

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

perform_phasing.sh 

Options:

-i <path of input directory where all chromosome folders are present>

-c <chr|CHR> Folder names start with either CHR or chr

-p <prefix of .bim .bed file> 

-r <path for reference files>

bash /home/ss5505/scripts_cumc/gwas_/phasing_scripts/perform_phasing.sh -i /home/ss5505/tool_toy_files/shapeit_impute2/test_data -r ~/reference_panels/1000GP_Phase3 -c CHR -p HGWAS6_final_

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

val=1239267
st="EGAZ0000"

for chr_num in {1..22}

do

    #--loop through chromosomes
    num=$(($val + $chr_num))
    printf "$st$num\n"
    
    #--reference files
    #--set bed file per chromsome and get hap, legend, map, sample files from reference panels
    b_file=$input_dir/$chrom_pre$chr_num/$prefix_input_files$chrom_pre$chr_num #store bed/bim/fam file
    if [ "$chr_num" -ne 1 ]
    then

	m_file=$ref_path/genetic_map_chr${chr_num}_combined_b37.txt
	hap_file=$ref_path/"$st$num"_HRC.r1-1.EGA.GRCh37.chr"$chr_num".hap.gz
	legend_file=$ref_path/"$st$num"_HRC.r1-1.EGA.GRCh37.chr"$chr_num".legend.gz
	sample_file=$ref_path/"$st$num"_HRC.r1-1.EGA.GRCh37.chr"$chr_num".samples
	
    else
	#CHR 1 for HRC is different
	#
	m_file=$ref_path/genetic_map_chr${chr_num}_combined_b37.txt
	hap_file=$ref_path/EGAZ00001239268_HRC.r1-1.EGA.GRCh37.chr1.noIBD.hap.gz
	legend_file=$ref_path/EGAZ00001239268_HRC.r1-1.EGA.GRCh37.chr1.noIBD.legend.gz
	sample_file=$ref_path/EGAZ00001239268_HRC.r1-1.EGA.GRCh37.chr1.noIBD.samples
	
    fi

    #----
    output_log=$input_dir/$chrom_pre$chr_num/$chrom_pre$chr_num".alignments" #output name for alignments logs

    #--submit job 
#bash /home/ss5505/scripts_cumc/gwas_/phasing_scripts/launch_job_chrm_ph.sh -b $b_file -m $m_file -f $hap_file -l $legend_file -s $sample_file -p $chrom_pre$chr_num

  qsub -wd "$input_dir/$chrom_pre$chr_num/" -l h_vmem=35G -N "ph_chr""$chr_num" -v bash /mnt/mfs/hgrcgrid/homes/ss5505/scripts_cumc/gwas_/phasing_scripts/HRC_phasing_impute2/launch_job_chrm_ph.sh -b $b_file -m $m_file -f $hap_file -l $legend_file -s $sample_file -p $chrom_pre$chr_num

done
