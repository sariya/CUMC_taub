#!/bin/bash

##############################
#Date 04/02/2019
#Sanjeev Sariya
##############################

for i in {12..21}

do

   grep -v "#" filtered_CHR"$i".vcf | cut -f 3,5  | grep "\."  | cut -f 1 > monomorphic_sites_CHR"$i"
   printf "CHR$i complete\n" 
done

for i in {1..9}

do

   grep -v "#" filtered_CHR"$i".vcf | cut -f 3,5  | grep "\."  | cut -f 1 > monomorphic_sites_CHR"$i"
   printf "CHR$i complete\n" 
done

for i in {12..21}

do

   val=($(grep -v "#" monomorphic_sites_CHR"$i" | cut -f 3,5  | wc -l))
   printf "CHR"$i" $val\n" >> STR_monomorphic_counts
done

for i in {1..9}

do
 val=($(grep -v "#" filtered_CHR"$i".vcf | cut -f 3,5  | wc -l))
   printf "CHR"$i" $val\n" >>  STR_counts 
done


for i in {1..9}

do

   val=($(grep -v "#" monomorphic_sites_CHR"$i" | cut -f 3,5  | wc -l))
   printf "CHR"$i" $val\n" >> STR_monomorphic_counts
done

