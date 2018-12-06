#!/bin/Rscript

#
#Sanjeev Sariya
#bin imputed SNPs into gene names
#
#Resize gene range based on your need
#

#
#have to work further when one snp falls into multiple gene boundaires
#

##Code from 
##https://jef.works/blog/2016/12/06/mapping-snps-and-peaks-to-genes-in-R/
## and biostars and support bio-conductor

##--gene file name
file.genes<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merged_PR_HGWAS/GWAS_analyses/output/gene_per_chr/genes_22.txt"
df.genes<-read.table(file.genes,header=TRUE)

row.names(df.genes)<-df.genes$Gene
df.genes<-df.genes[,-c(4)]

range2GRanges <- function(df) {
    require(GenomicRanges)
    require(IRanges)
	gr <- GenomicRanges::GRanges(
        seqnames = df[,1],
        ranges=IRanges(start = df[,2], end = df[,3])
        )
    return(gr)
}


string2range <- function(pos, delim=' ', region=TRUE) {
    posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
    posp[,1] <- posp[,1]
	posp[,2] <- as.numeric(as.character(posp[,2]))
	if(region) {
        posp[,3] <- as.numeric(as.character(posp[,3]))
	} else {
	    posp[,3] <- posp[,2]
	}
    return(posp)
}

gene.granges<-range2GRanges(df.genes)
names(gene.granges)<-row.names(df.genes)

#
#Extend ranges
#https://stackoverflow.com/questions/34331485/extend-range-in-both-directions
#

start(gene.granges)<-start(gene.granges) - 10000
end(gene.granges)<-end(gene.granges) + 10000


snps <- c("22:51164259:A:G", "22:51185517:A:G", "22:51232086:A:G")
snps.ranges <- string2range(snps, delim=":", region=FALSE)
snps.ranges<-snps.ranges[,-(4)]


snps.granges <- range2GRanges(snps.ranges)
names(snps.granges) <- snps

r1 <- snps.granges
r2 <- gene.granges


overlap <- GenomicRanges::findOverlaps(r1, r2)
# make vector of SNPs to genes

hits<-names(r2)[subjectHits(overlap)]

names(hits)<-names(r1)[queryHits(overlap)]

#hits <- names(r2)[slot(overlap, "subjectHits")]
#names(hits) <- names(r1)[slot(overlap, "queryHits")]
#hits


r1[names(hits),]

r2[hits,]

for (i in 1:length(hits)){
print(paste(names(hits[i]),hits[i]))

}

"
sessionInfo()
R version 3.4.2 (2017-09-28)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 9 (stretch)

Matrix products: default
BLAS: /mnt/mfs/cluster/bin/R-3.4/lib/libRblas.so
LAPACK: /mnt/mfs/cluster/bin/R-3.4/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
[1] GenomicRanges_1.30.3 GenomeInfoDb_1.14.0  IRanges_2.12.0
[4] S4Vectors_0.16.0     BiocGenerics_0.24.0

loaded via a namespace (and not attached):
[1] zlibbioc_1.24.0        compiler_3.4.2         tools_3.4.2
[4] XVector_0.18.0         GenomeInfoDbData_1.0.0 RCurl_1.95-4.11
[7] bitops_1.0-6

"

