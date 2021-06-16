#Ssanjeev Sariya
#06 16 2021



CairoTIFF(filename = "coloredbox_brains.tiff", height=800, width=900) ## TIFF image
##CairoTIFF(filename = "box_brains_04.26.2021.tiff", height=80, width=80, units="mm", dpi=300) ## TIFF image

ggplot(prs_df_chbrains, aes(x=status,y=RNOmni::RankNorm(final_sum, k=3/8), color=status)) + 
geom_boxplot(lwd=0.45, width =0.25) + theme_minimal() + 
labs(title="", x="Status", y="Polygenic Risk Score") + 
theme( plot.title = element_text( size=10, face="bold", hjust = 0.5),
axis.title.x = element_text(size=7, face="bold"),
axis.title.y = element_text(size=7, face="bold"),
axis.text.y = element_text(size=6,face="bold"),
axis.text.x=element_text(size=6,face="bold"),
legend.key.size = unit(2.0, "mm"),
legend.background = element_rect(size=0.3,linetype="solid", colour ="black"), 
legend.text=element_text(size=4), 
legend.title=element_text(size=5),
legend.position = c(0.9, 0.8) ) + 
labs(color='Status') 

dev.off()

##############################


prs_df_chbrains<-read.table("CHbrains_N33_phenotype_PRS.txt",header=TRUE)
CairoTIFF(filename = "gray_colored.tiff", height=80, width=80, units="mm", dpi=300) ## TIFF image

ggplot(prs_df_chbrains, aes(color=status, x=status,y=RNOmni::RankNorm(final_sum, k=3/8) )) + 
geom_boxplot(lwd=0.45, width =0.25) + scale_color_grey()+
theme_minimal() +  
labs(title="",x="Status",y="Polygenic Risk Score") + 
theme( plot.title = element_text(size=10, face="bold", hjust = 0.5),
axis.title.x = element_text(size=7,face="bold"),
axis.title.y = element_text(size=7,face="bold"),
axis.text.y = element_text(size=6,face="bold"),
axis.text.x=element_text(size=6,face="bold"),
legend.key.size = unit(2.0,"mm"),
legend.background = element_rect(size=0.3,linetype="solid", colour ="black"), 
legend.text=element_text(size=4), 
legend.title=element_text(size=5),
legend.position = c(0.9, 0.8) ) + 
labs(color='Status') 

dev.off()

