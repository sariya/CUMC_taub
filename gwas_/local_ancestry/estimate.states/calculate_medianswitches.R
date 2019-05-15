#!/usr/bin/Rscript

#
#05/15/2019
#Sanjeev Sariya
#use this script after output from recode_LA_haplotypes.R
#memory request was 35G, used 9.X G.

#input individuals are 10,758
#I got this code from Daniel Shriner 

library(coda)

res <- matrix(data=NA,nrow=10758,ncol=22) ##store value per person per CHR

for (n in 1:22) {
    
    dat <- read.table(paste("recode_CHR",n,"_states.txt",sep=""))

    for (j in 1:10758) {
        res[j,n] <- effectiveSize(na.omit(unlist(dat[,j])))
    }
    ##for loop ends for per person
    
    print(paste("Done with CHR",n,sep=""))
}
##for loop ends
 
median(apply(res,1,sum))
##[1] 318.787

print("Complete")

##R version 3.4.2 (2017-09-28)
##coda_0.19-2
