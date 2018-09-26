#!/bin/bash

#Sanjeev Sariya
#09/26/2018
#Identify files with varying number of columns
#Loop through output of GMMAT chunk. Print file names that have different columns. 
#

chunk_files=($(ls CHR*_chunk*_scoremodel2))

for x in "${chunk_files[@]}"
do
#printf "$x\n"
chunk_name=($(basename $x) ) 

count=($(awk '{print NF}' $x | sort | uniq | wc -l ))
#printf "number of lines are $count\n"

if [ "$count" -ne 1 ] 
then 
printf "$x has different number of columns\n"
fi

done
