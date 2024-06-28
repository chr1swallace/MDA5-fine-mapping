
## * reference data
## get 1000 Genomes ref data, both the original, and the updated "30X" version.  They have slightly different SNP coverage.

## get ukbb data - NB build 37
BBDIR="/home/cew54/biobank/genotypes-imputed"
input=paste0(BBDIR,"/ukb22828_c2_b0_v3.bgen")
if(!file.exists("ukbb_hg37.bim")) {
    cmd=paste0("bgenix -g ",input," -list > snps_hg37.csv && ",
               "awk \'$0 !~ /#/ {print \"chr\"$3, $4, $4+1, $1}\' snps_hg37.csv | sed '1d' | sed 's/chr0/chr/g' > ~/scratch/tmp_c2.bed && ",
               "liftOver ~/scratch/tmp_c#{chr}.bed /home/cew54/share/Data/reference/hg19ToHg38.over.chain.gz ~/scratch/tmpout_c2.bed ~/scratch/tmpunmap_c2.bed\n")
    cat(cmd,file="torun.sh")
    system("qlines.rb -r torun.sh") # runs the above command on HPC slurm queue
    message("respond when above job is complete")
    a=readLines(n=1) 
    system("./ukbb_merge_snpbuilds.R")
    system("rm ~/scratch/tmp_2.bed ~/scratch/tmpout_2.bed ~/scratch/tmpunmap_2.bed")
}

## read in snps for alignment of gwas data
snps=fread("ukbb_hg38.bim")
setnames(snps, c("chr","rsid","null","position","a0","a1"))
snps[,snp:=paste(chr,position,a0,a1,sep="-")]
snps[,asnp:=paste(chr,position,a1,a0,sep="-")]
snps37=fread("ukbb_hg37.bim")
setnames(snps37, c("chr","rsid","null","position","a0","a1"))
snps[,pos37:=snps37$position]

## * load GWAS summary data and run SuSiE for fine mapping - can skip this section if already run

if(!file.exists("gwas.RData")) {
    source("run-gwas.R")
} else {
    (load("gwas.RData"))
}


variants=lapply(DATA, "[[", "snp")  %>% unlist()  %>% unique()

## run extract_....sh here, filtering on variants (rsid) from above before doing plink
if(!file.exists("ukbb40k_hg37.bim")) {
    hg37=fread("ukbb_hg37.bim")
    hg38=fread("ukbb_hg38.bim")
    setnames(hg37,c("chr","rsid","null","position","a0","a1"))
    setnames(hg38,c("chr","rsid","null","position","a0","a1"))
    hg37[,snp:=paste(chr,position,a0,a1,sep="-")]
    hg38[,snp:=paste(chr,position,a0,a1,sep="-")]
    m=merge(hg37,hg38,by="rsid",suffixes=c(".37",".38"))
    m=m[a0.37==a0.38 & a1.37==a1.38]
    rsid=m[snp.38 %in% variants]$rsid
    cat(rsid,file="rsid.tmp",sep="\n")
    system("./extract_ukbb40k.sh") # creates ukbb40k_hg37.bim containing 40k eur from ukbb
    ## cat("./extract_ukbb40k.sh\n", file="torun.sh") # creates ukbb40k_hg37.bim containing 40k eur from ukbb
    ## system("qlines.rb -r torun.sh") # runs the above command on HPC slurm queue
}

 
if(!file.exists("ld-maf.RData")) {
    stub="ukbb40k_hg37"
    paste0("plink --r square --keep-allele-order --bfile ",stub," --out ",stub)  %>% system()
    LD=paste0(stub,".ld")  %>% scan(., what=0)
    n=length(LD)  %>% sqrt()
    LD  %<>%  matrix(., n, n, byrow=TRUE)
    paste0(stub,".ld")  %>% unlink()
    paste0("plink --freq --keep-allele-order --bfile ",stub," --out ",stub)  %>% system()
    MAF=paste0(stub,".frq")  %>% fread()
    paste0(stub,".frq")  %>% unlink()
    MAFsnps=paste0(stub,".bim")  %>% fread()
    MAF=MAF$MAF
    names(MAF)=with(MAFsnps, paste0(V1,"-",V4,"-",V5,"-",V6))
    ## move to hg38
    if(!exists("m")) {
       hg37=fread("ukbb_hg37.bim")
       hg38=fread("ukbb_hg38.bim")
       setnames(hg37,c("chr","rsid","null","position","a0","a1"))
       setnames(hg38,c("chr","rsid","null","position","a0","a1"))
       hg37[,snp:=paste(chr,position,a0,a1,sep="-")]
       hg38[,snp:=paste(chr,position,a0,a1,sep="-")]
       m=merge(hg37,hg38,by="rsid",suffixes=c(".37",".38"))
    } 
    mm=match(names(MAF),m$snp.37)
    names(MAF)=m$snp.38[mm]
    dimnames(LD)=list(names(MAF),names(MAF))
    LD[is.na(LD)]=0
    use=MAF<1 & MAF > 0
    MAF=MAF[use]
    LD=LD[use,use]
    save(LD,MAF,file="ld-maf.RData")
} else {
    (load("ld-maf.RData"))
}
