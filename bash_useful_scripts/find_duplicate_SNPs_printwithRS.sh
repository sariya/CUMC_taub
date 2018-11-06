#!/bin/bash

#
#Sanjeev Sariya
#11/05/2018

#--if first column starts with rs 
awk '{

if( ($2 ~ /^rs/)){ print $0
} }' RS_IDS_CHR9 > only_RSidsCHR9

###-----------------------------

#----print duplicates $1- is chr:pos SNP id

awk ' {
if($1 in ARRAY){
ARRAY[$1]=ARRAY[$1]+1
}
else{
ARRAY[$1]=1
}
}

END{
for (key in ARRAY){

if(ARRAY[key] >1){
#--if value more than 1 print

print key,ARRAY[key]
} } } ' RS_IDS_CHR9 

