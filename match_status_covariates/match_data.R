#!/usr/bin/env Rscript

#
#Date 01/02/2020
#Sanjeev Sariya
#

check_columns<-function(cols1,cols2){

print("we are here")
if(length(cols1)!=length(cols2)){
stop("column number isn't same. Exiting")

}
##if ends to check column number

if(length(match(colnames(df2),colnames(df1))) != length(colnames(df2))){
stop("columns don't match. Exiting")

}
##first match ends 

if(length(match(colnames(df1),colnames(df2))) != length(colnames(df1))){
stop("columns don't match. Exiting")

}
##second match ends 


}

#function check_columns ends 

file1<-"input_one.txt"
file2<-"input_two.txt"

df1<-read.table(file1,header=TRUE)
df2<-read.table(file2,header=TRUE)

check_columns(colnames(df1),colnames(df2))


