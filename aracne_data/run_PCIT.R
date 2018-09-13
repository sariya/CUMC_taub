#!/bin/Rscript

#
#
#Date 09 1 2018
#Running PCIT on aracne data
#
library(Cairo)
library("PCIT")
setwd("/mnt/mfs/scratch/GT_Admix/ss5505/gene_networks/PCIT/genes_mapped")

df_expression<-data.table::fread("/mnt/mfs/scratch/GT_Admix/ss5505/gene_networks/PCIT/genes_mapped/complete_genes_mapped")
print(dim(df_expression))
print(df_expression[1:10,1:10])
df_expression_transformed<-t(df_expression)
print("Transformed")

colnames(df_expression_transformed)<-df_expression_transformed[1,]
df_expression_transformed<-df_expression_transformed[-c(1),]
df_expression_transformed<-apply(df_expression_transformed, 2, as.numeric) #make values numeric

correlation_df<-cor(df_expression_transformed)
print(correlation_df[1:10,1:10])

print("Correlation found")

# write correlation table on file. Use it later for plots

write.table(correlation_df,"correlation_matrix.txt",sep = "\t",quote=FALSE,col.names = TRUE,row.names = FALSE)

system.time(result_serial <- pcit(correlation_df, force.serial=TRUE))
system.time(result <- pcit(correlation_df))

##
# pcit() doesn't return the ind's in the same order when done in parallel and serial
# check that we got the same answer using both functions

all.equal(correlation_df[idx(result_serial)], correlation_df[idx(result)])

# get the matric indices for the meaningful and unmeaningful correlations
meaningful.idx <- idx(result_serial)

unmeaningful.idx <- idxInvert(nrow(correlation_df), meaningful.idx)

# create a copy of the correlation matrix and set unmeaingful correlations to zero
correlation_df.new <- correlation_df

correlation_df.new[unmeaningful.idx] <- 0

# convert adjacency matrix into an edge list
edgeList <- getEdgeList(correlation_df.new)

#--print edgelist values in a file. use for cytoscape
write.table(edgeList ,"edge_list_weight",sep = "\t",quote=FALSE,col.names = TRUE,row.names = FALSE)

cc <- clusteringCoefficient(correlation_df.new)
ccp <- clusteringCoefficientPercent(correlation_df.new)

print("CC percent is ")
print(ccp)
write.table(cc ,"clustering_coeff.txt",sep = "\t",quote=FALSE,col.names = TRUE,row.names = FALSE)

#--make 
CairoJPEG("all_plots_matrix.jpeg",quality=75,height=900,width=1200)
op <- par(mfrow=c(3,2))

plot(density(correlation_df[upper.tri(correlation_df)]), main="Density Plot of Raw Correlation Coefficients", xlab="Correlation Coefficient")

hist(cc, main="Connectivity Distribution", xlab="Proportion of Connections", ylab="Number of Genes")
hist(cc*length(cc), main="Connectivity Distribution", xlab="Number of Connections", ylab="Number of Genes")

# plot the distribution of all correlations superimposed by that of the meaningful corrections in black
plotCorCoeff(correlation_df, list("PCIT Significant" = meaningful.idx), col=c("black"))

# plot the distribution of all correlations superimposed by that of the meaningful connections in black and the absolute correlations > 0.5 in red
abs.idx <- which(abs(correlation_df)>  0.5)

# we'll change the order and use some transparent colours using rgb()
plotCorCoeff(correlation_df, list("PCIT Significant" = meaningful.idx, "abs. cor. > 0.5" = abs.idx), col=c(rgb(1,0,0,0.7), rgb(0,0,0,0.7)))

par(op)
dev.off()

print("Done iwth plotting")
print("Exiting!!")

"sessionInfo()
R version 3.4.2 (2017-09-28)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 9 (stretch)

Matrix products: default
BLAS: /mnt/mfs/cluster/bin/R-3.4/lib/libRblas.so
LAPACK: /mnt/mfs/cluster/bin/R-3.4/lib/libRlapack.so

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] Cairo_1.5-9 PCIT_1.5-3

loaded via a namespace (and not attached):
[1] compiler_3.4.2    tools_3.4.2       data.table_1.11.4"
