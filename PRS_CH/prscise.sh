#!/usr/bin/bash

PRSice_linux  --base  GWAS_filteredAPOE \
--target CHR#_PRS--nonfounder \
--bar-levels 0.00000005,0.00001,0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5  \
--no-full --fastscore \
--pheno phenotype_PRS.txt \
--pheno-col AD \
--binary-target T \
--cov phenotype.txt \
--cov-col PC1,PC2,PC3,PC4 \
--print-snp --seed 2397373689 \
--thread 8 --pvalue PVAL --snp SNP --bp BP --chr CHR --stat OR \
--A2 Allele1 --A1 Allele2 --model add --type bed --score avg \
--clump-kb 250 --clump-r2 0.1 --clump-p 1 --maf 0.01 \
--perm 10000 --logit-perm --all-score --out PRS_final

# PRScise version 2.3.3 (2020-08-05)
