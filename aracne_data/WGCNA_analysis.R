#!/bin/Rscript

#
#update 02 17 2018
#Sanjeev Sariya
#
#Perform WGCNA using data from Dr. Amanda Myers. 
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

#
#We work with Amanda myers data 363 individuals. We clean samples and we clean genes
#
#Tutorial at: #http://pklab.med.harvard.edu/scw2014/WGCNA.html##################
#

library(dplyr)
library(WGCNA)
library(Cairo)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

exprs_data<-read.table("duplicate_genes_summed",header=TRUE)
print(dim(exprs_data))

#transpose the expression data
data_exprs.cleaned<-as.data.frame(t(exprs_data[, -c(1)])); #remove gene column

#add row names, and col names
names(data_exprs.cleaned) = exprs_data$gene
rownames(data_exprs.cleaned) = names(exprs_data)[-c(1)]
print(dim(data_exprs.cleaned))

print(data_exprs.cleaned[1:10,1:10])

#check data for excessive missing values and identi_cation of outlier microarray
gsg = goodSamplesGenes(data_exprs.cleaned, verbose = 3);

#--everything OK with mapped genes
if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(data_exprs.cleaned)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(data_exprs.cleaned)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
data_exprs.cleaned= data_exprs.cleaned[gsg$goodSamples, gsg$goodGenes]
}

#Check outliers
sampleTree = hclust(dist(data_exprs.cleaned ), method = "average");

# Plot the sample tree: 
# The user should change the dimensions if the window is too large or too small.

#Not removing any samples from Amanda myers brain expression data
#
CairoJPEG("sample120_outliers_tree.jpeg",width=1200,height=900)
par(cex = 0.8);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
abline(h=120, col = "red")
dev.off()

print("Printed JPEG file for sample tree to see outlier")

#-for amanda myers data. I'm using 140 as cut height and losing 29 samples out of 363. Less than 10% is acceptable
labels_def = cutreeStatic(sampleTree, cutHeight = 140)  #--default cluster size is 50
table(labels_def )
keepSamples = (labels_def==1)

datExpr = data_exprs.cleaned [keepSamples, ]
nGenes = ncol(data_exprs.cleaned )
nSamples = nrow(data_exprs.cleaned )

#datExpr  is new data expression frame
##################################################################
#Automatic, one-step network construction and module detection:
##################################################################

#datExpr  is new data expression frame
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) # work with rmoeved samples
CairoJPEG("soft_gene.jpeg",height=900,width=1000)

par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()
#############################

softPower = 6; #for amnda myers data we are using 6 as sft threshold power 334 samples
adjacency = adjacency(datExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

print(dim(dissTOM)) # -- should be sq matrix with genes

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
CairoJPEG("result_bnwclustering_gene_0.4.jpeg",height=900,width=1000)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)

#
#21422 genes
#
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

CairoJPEG("module_colors.jpeg",height=900,width=1000)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

#discard the unassigned genes, and focus on the rest
restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = softPower) #samples cleaned

print(dim(diss1)) #19499 19499 we lose 1923 genes
colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
hier1=hclust(as.dist(diss1), method="average" )
CairoJPEG("module_colors_nogrey.jpeg",height=900,width=1000)
plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()


#set the diagonal of the dissimilarity to NA 
diag(diss1) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
CairoJPEG("module_colors_nogrey_heatmap.jpeg",height=900,width=1000)
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]),main="Module clustered")
dev.off()

###Extract modules

length(which(dynamicMods == 0)) #0 means grey color 
datExpr_cleaned_greycolor<-datExpr[,-c(which(dynamicMods == 0))]
gene.names<-colnames(datExpr_cleaned_greycolor) 
#gene.names<-colnames(datExpr)
SubGeneNames<-gene.names


dynamiccolors_nogrey<-dynamicColors[dynamicColors!="grey"]

module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
    module=SubGeneNames[which(dynamiccolors_nogrey==color)]
    write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

module.order <- unlist(tapply(1:ncol(datExpr_cleaned_greycolor),as.factor(dynamiccolors_nogrey),I)) #samples cleaned

m<-t(t(datExpr_cleaned_greycolor[,module.order])/apply(datExpr_cleaned_greycolor[,module.order],2,max)) #samples cleaned

CairoJPEG("heatmap_expression_profiles.jpeg",height=900,width=1000)
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])
dev.off()

### get value per moduel for person
###Quantify module similarity by eigengene correlation. Eigengenes: Module representatives
####

MEList = moduleEigengenes(datExpr_cleaned_greycolor, colors = dynamiccolors_nogrey)
##MEList = moduleEigengenes(data_exprs.cleaned, colors = dynamicColors)
MEs = MEList$eigengenes

CairoJPEG("eigennetworks.jpeg",height=900,width=1000)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
dev.off()

##############################################################################################

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
CairoJPEG("clustering_module_eigengenes_cutoff_30.jpeg",height=900,width=1000)

plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=0.30,col="red") #plot line to merge similar modules
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(datExpr_cleaned_greycolor, dynamiccolors_nogrey, cutHeight = 0.30, verbose = 3) #we have 40 gene modules
length(merge$colors) #19499 - after excluding grey color

##merge = mergeCloseModules(data_exprs.cleaned, dynamicColors, cutHeight = 0.30, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#read disease status
diseasestatus<-read.table("disease_status_sample",header=TRUE)
sample_cleaned_expression<-rownames(datExpr_cleaned_greycolor) ## get sample names

#merge them
merge_samples_diseasestatus<-(merge(rownames(datExpr_cleaned_greycolor),diseasestatus,by=0)) 

#make row.names as numeric 
merge_samples_diseasestatus$Row.names <- as.numeric(merge_samples_diseasestatus$Row.names)

#sort based on 
merge_samples_diseasestatus<-(merge_samples_diseasestatus[order(merge_samples_diseasestatus$Row.names),] )
 

###Run logistic regression to check significance of any module

df_colors_pvalue <- data.frame(matrix(ncol = 2, nrow =length(colnames(mergedMEs) )))
colnames(df_colors_pvalue)<-c("colorModule","Pvalue")


#
#mergedMEs it has column count as colors and row count as samples 
#
for (i in 1:length(colnames(mergedMEs))){
##=- loop over different module's eigen values

#--run linear logistic regression
temp.pvalue<-(summary(lm(formula = merge_samples_diseasestatus$Status ~ mergedMEs[,i], data = MEs, family = "binomial"))$coeff)[2,4]
df_colors_pvalue[i,1] <-colnames(mergedMEs)[i]
df_colors_pvalue[i,2] <-temp.pvalue

}

df_colors_pvalue<- df_colors_pvalue[order(df_colors_pvalue$Pvalue),] 


write.table(df_colors_pvalue,"colors_pvalue_regression",sep = "\t",quote=FALSE,col.names = TRUE,row.names = FALSE)

#--use manual's method
significance_module<-t(as.data.frame(signif(cor(merge_samples_diseasestatus$Status,mergedMEs, use="p"),2)))

significance_module<- significance_module[order(significance_module[,1]),] 
significance_module<-as.data.frame(significance_module)
significance_module$colors<-row.names(significance_module)
write.table(significance_module,"colors_significance",sep = "\t",quote=FALSE,col.names = TRUE,row.names = FALSE)


##--just verify things again
p.values = corPvalueStudent(cor(diseasestatus$Status,mergedMEs, use="p"), nSamples = length(diseasestatus$Status))

# Measure of module significance as average gene significance

GS1=as.numeric(cor(diseasestatus$Status,data_exprs.cleaned, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.

ModuleSignificance=tapply(GeneSignificance, merge$colors, mean, na.rm=T)
ModuleSignificance<-(as.data.frame(ModuleSignificance))
ModuleSignificance$color<-row.names(ModuleSignificance)

ModuleSignificance<- ModuleSignificance[,c(2,1)]

CairoJPEG("merged_colors_threshold.jpeg",height=900,width=1000,quality=75)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

#In the subsequent analysis, we will use the merged module colors in mergedColors. We save the relevant variables for
#use in subsequent parts of the tutorial:

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

##################################################################
## one step network construction
net = blockwiseModules(data_exprs.cleaned, 
power = 6,TOMType = "unsigned", minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,saveTOMFileBase = "Sample_genes_mapped",verbose = 3)

##############################################################################################
kIM = intramodularConnectivity(adjacency, dynamicColors, scaleByMax = TRUE) 

