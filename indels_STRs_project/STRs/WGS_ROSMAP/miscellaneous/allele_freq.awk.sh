#!/bin/bash

#
#Date 04/20/2019
#Sanjeev Sariya
##Local script to get allele frequency of grepped STRs from in-house dosage created.

 awk -F"," '{
for (i=4;i<=NF;i++){

if($i!="NA"){
total=total+1
if($i!=0){
sum=sum+$i

}}}

print sum/(total*2) 
}'