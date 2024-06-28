#!/usr/bin/env Rscript

library(randomFunctions)
library(data.table)
library(magrittr)

infile="snps_hg37.csv"
mapped=paste0("~/scratch/tmpout_c2.bed")

st=161267074
en=163318684

snps=fread(infile)
hg38=fread(mapped)
m=match(hg38$V4,snps$alternate_ids)
snps$hg38=as.numeric(NA)
snps$hg38[m[!is.na(m)]] = hg38$V2[!is.na(m)]
snps=snps[ !is.na(hg38) ] # & number_of_alleles==2 ]
snps=snps[ hg38>=st & hg38<=en ]
## snps=snps[ !is.na(hg38) & minor_allele_frequency > 0.005 & HW_lrt_p_value > 1e-6 & is.finite(HW_lrt_p_value) &
##           nchar(alleleA)==1 & nchar(alleleB)==1]

fwrite(snps, file="snps_hg38.csv")
## also create plink bims for ifih1 block

snps$zero=0
fwrite(snps[ , .(chromosome, rsid, zero, hg38, first_allele, alternative_alleles)],
       file="ukbb_hg38.bim",
       col.names=FALSE, sep="\t")
fwrite(snps[ , .(chromosome, rsid, zero, position, first_allele, alternative_alleles)],
       file="ukbb_hg37.bim",
       col.names=FALSE, sep="\t")
