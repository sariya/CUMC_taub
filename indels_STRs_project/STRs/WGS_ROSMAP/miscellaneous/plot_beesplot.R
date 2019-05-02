#!/bin/RScript

#
#Date 04/24/2019
#

library(ggbeeswarm)
library(Cairo)
file.input<-"globalpathology.input" ##this is case vs controls 
df.input<-read.table(file.input,header=TRUE)
print(dim(df.input))


CairoJPEG(filename = "globalPath_AD_nonAD.jpeg",quality = 75,width=1200,height=900) 

ggplot(df.input,aes(pathoAD,gpath,color=factor(pathoAD))) + geom_beeswarm() +
ggtitle("Global pathology for Pathological AD vs nonAD ") + labs(x = "Pathological AD",y="Global Pathology") +
labs(caption = "(based on 1,134 individuals)") +theme(plot.title = element_text(size=20, face="bold.italic"),
axis.title.x = element_text(size=18, face="bold"),
axis.title.y = element_text( size=18, face="bold"),
axis.text=element_text(size=16),legend.text=element_text(size=14),
plot.caption = element_text(color = "blue", face = "italic",size=12)) 


CairoJPEG(filename = "globalPath_AD_nonAD.jpeg",quality = 75,width=1200,height=900) 

ggplot(df.input,aes(pathoAD,gpath,group=pathoAD,colour=pathoAD)) + geom_beeswarm() +
ggtitle("Global pathology for Pathological AD vs nonAD ") + labs(x = "Pathological AD",y="Global Pathology") +
labs(caption = "(based on 1,134 individuals)") +theme(plot.title = element_text(size=20, face="bold.italic"),
axis.title.x = element_text(size=18, face="bold"),
axis.title.y = element_text( size=18, face="bold"),
axis.text=element_text(size=16),legend.text=element_text(size=14),
plot.caption = element_text(color = "blue", face = "italic",size=12))  +  guides(color=guide_legend("Status"))

dev.off()



