#!/bin/bash

#
#Date 04/20/2019
#Sanjeev Sariya
##Local script to get allel frequency of grepped STRs from in-house dosage created.

awk -F"," '{
for (i=4;i<=NF;i++){
if($i!="NA"){
total=total+1
if($i!=0){
minor=minor+1

}
}

}

print minor,total,minor/(2*total)
}'
