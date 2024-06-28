#!/bin/bash
BBDIR="/home/cew54/biobank/genotypes-imputed"
#chr=`echo $block | sed 's/chr//; s/_block[0-9]*//; s/a//'`
chr='2'

bgen_samples="${BBDIR}/ukb22828_c${chr}_b0_v3_s487160.sample"
bgen_file="${BBDIR}/ukb22828_c${chr}_b0_v3.bgen"

## list of samples
samples="EUR40k.sample"
## list of snps
tmp_snps="rsid.tmp"

## bgenix to extract the snp subset, then qctool to extract the samples and convert to plink format
## tmp_bgen=tempfile(fileext=".bgen")
tmp_bed=`mktemp`
bgenix -g ${bgen_file} -incl-rsids ${tmp_snps} | \
    qctool_v2.2.0 -g - -filetype bgen -s ${bgen_samples} -incl-samples ${samples} \
		  -og ${tmp_bed} -ofiletype binary_ped  

## plink to clean up
out_bed="ukbb40k_hg37" ## qctool will also make .bim and .fam
plink --geno 0.01 --snps-only --hwe 0.00000001 --make-bed \
      --keep-allele-order --bfile ${tmp_bed} --out ${out_bed}
rm ${tmp_bed}.*

## plink with alt_ids
out_bim="ukbb40k_hg37.bim" ## qctool will also make .bim and .fam
alt_bim="ukbb40k_alt_hg37.bim" ## qctool will also make .bim and .fam
cat ${out_bim} | awk 'BEGIN{ OFS="\t" } {print $1,$1 ":" $4 "_" $5 "_" $6,$3,$4,$5,$6}' > ${alt_bim}
	 
