#!/bin/Rscript

#
#
#Date 09 1 2018
#Running PCIT on aracne data
#
library("PCIT")

df_expression<-data.table::fread("duplicate_genes_summed")
print(dim(df_expression))
print(df_expression[1:10,1:10])
df_expression_transformed<-t(df_expression)
print("Transformed")

gene_names<-df_expression_transformed[1,]
df_expression_transformed<-df_expression_transformed[-c(1),]


df_expression_transformed<-apply(df_expression_transformed, 2, as.numeric)

correlation_df<-cor(df_expression_transformed)
print(correlation_df[1:10,1:10])

print("Correlation found")

system.time(result_serial <- pcit(correlation_df, force.serial=TRUE))





"sessionInfo()
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
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] PCIT_1.5-3

loaded via a namespace (and not attached):
[1] compiler_3.4.2"