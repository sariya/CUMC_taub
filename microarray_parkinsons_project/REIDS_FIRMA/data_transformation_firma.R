#!/bin/Rscript

#
#
#Sanjeev Sariya
#11/14/2018
#Gius - KM PD expression data

#
#Read output from FIRMA REIDS and transform to log
#

load("CEL_FIRMA_Output.RData")
print(dim(CEL_FIRMA_Output))

log_transformed_firmavalues<-log(CEL_FIRMA_Output[,3:ncol(CEL_FIRMA_Output)])
temp_exoncolnames<-colnames(log_transformed_firmavalues)

#--fix col names

for(i in 1:length(temp_exoncolnames)){

temp_exoncolnames[i]<-(paste("X",temp_exoncolnames[i],sep=""))
temp_exoncolnames[i]<-(paste(temp_exoncolnames[i],".CEL",sep=""))
temp_exoncolnames[i]<-gsub("-",".",temp_exoncolnames[i])

}

#################################
#bind columns
#################################

colnames(log_transformed_firmavalues)<-temp_exoncolnames
bind_exongeneids_log<- cbind(CEL_FIRMA_Output[,1:2],log_transformed_firmavalues)
list_split_firmaoutput<-split(bind_exongeneids_log, (seq(nrow(bind_exongeneids_log))-1) %/% 25000) 
print(length(list_split_firmaoutput))

write.table(bind_exongeneids_log, "firma.exons.txt", sep="\t", quote=F, row.names=FALSE,col.names=TRUE)

for(i in 1:length(list_split_firmaoutput) )
{
file<-paste("firma_exonchunks",i,sep="_") #--make file name for chunk
file_dest<-paste("./",file,sep="/")
write.table(list_split_firmaoutput[[i]], file_dest, sep="\t", quote=F, row.names=FALSE,col.names=TRUE)
}

print(dim(bind_exongeneids_log)) #284258    226

