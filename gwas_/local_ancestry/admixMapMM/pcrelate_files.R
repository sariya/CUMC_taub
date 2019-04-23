#!/usr/bin/Rscript

#
#Date 04/23/2019
#Sanjeev Sariya
#Local Ancestry
#https://rdrr.io/bioc/GENESIS/man/admixMapMM.html
#
#We need pcrelate data for admixmapMM. Use this script for .gds file creation in intermediatery steps
#
print("Begin Script")

library("SNPRelate")
library("GENESIS")
library("GWASTools")

print("loaded libraries")
