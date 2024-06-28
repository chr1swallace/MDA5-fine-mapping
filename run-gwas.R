 toget=c(liu_uc="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003045/harmonised/26192919-GCST003045-EFO_0000729.h.tsv.gz",
            liu_cro="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003044/harmonised/26192919-GCST003044-EFO_0000384.h.tsv.gz",
            lange_uc="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004133/harmonised/28067908-GCST004133-EFO_0000729.h.tsv.gz",
            lange_cro="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004132/harmonised/28067908-GCST004132-EFO_0000384.h.tsv.gz",
            forgetta_t1d="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010681/harmonised/32005708-GCST010681-EFO_0001359.h.tsv.gz",
            chiou_t1d="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/harmonised/34012112-GCST90014023-EFO_0001359-Build38.f.tsv.gz",
            robs_t1d_eur="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013445/GCST90013445_buildGRCh38.tsv",
            robs_t1d="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013446/GCST90013446_buildGRCh38.tsv",
            tsoi_pso="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005527/harmonised/23143594-GCST005527-EFO_0000676.h.tsv.gz")
            ## sakaue_hypo="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018862/harmonised/34594039-GCST90018862-EFO_0004705.h.tsv.gz"

    if(any(duplicated(toget)))
        stop("toget contains dups")
    
    d="~/scratch/ifih1"
    DATA=vector("list",length(toget))
    names(DATA)=names(toget)
    for(i in seq_along(toget)) {
        cat(i,"\t")
        target=toget[i]
        local=file.path(d,basename(target))
        if(!file.exists(local))
            system(paste("cd ",d," && wget ",target))
        tmp=fread(local)
        if("hm_beta" %in% names(tmp)) { # harmonised
            tmp[,pos38:=hm_pos][,chromosome:=hm_chrom][,P:=p_value][,REF:=hm_other_allele][,ALT:=hm_effect_allele]
            ## can't do this, because we end up with columns with duplicate names
            ## setnames(tmp,
            ##          c("hm_pos","hm_chrom","p_value","hm_other_allele","hm_effect_allele"),
            ##          c("pos38","chromosome","P","REF","ALT"))
            tmp[,beta:=hm_beta]
        } else {
            setnames(tmp,
                     c("base_pair_location","p_value","other_allele","effect_allele"),
                     c("pos38","P","REF","ALT"))
        }
        tmp[,disease:=sub("[a-z]*_","",names(toget)[i])  %>% sub("_eur","",.)]
        print(names(tmp))
        DATA[[i]]=tmp[chromosome==2 & pos38 >= st & pos38 <= en][,study:=names(toget)[i]]
        rm(tmp)
    }
    save(DATA, file="gwas.RData")

## add hypothyroidism from UKBB
f="~/rds/rds-phewas-kJB6O1vczqQ/D4_PanUKBB/summ_stats/categorical-20002-both_sexes-1226.tsv.bgz"
paste0("cp ",f," ~/scratch/ifih1/ukbb_hypo.tsv.gz")  %>% system()
d=fread("zcat ~/scratch/ifih1/ukbb_hypo.tsv.gz")
head(d)
d=d[chr==2]
d=merge(d,snps[,.(pos=pos37,position)],by="pos")

DATA$ukbb_hypo=d[,.(chromosome=chr,pos38=position,REF=ref,ALT=alt,beta=beta_EUR,standard_error=se_EUR,
                    P=10^(-neglog10_pval_EUR),study="ukbb_hypo",disease="hypo",
                    effect_allele_frequency=af_controls_EUR,
                    snp=paste(chr,pos,ref,alt,sep="-"))]


## LD and MAF
## get all variants to filter first
DATA=lapply(DATA, function(d) {
    tmp=d[chromosome==2 & pos38 > st & pos38 < en & !is.na(beta)]
    tmp[,snp:=toupper(paste(chromosome,pos38,REF,ALT,sep="-"))]
    tmp1=tmp[snp %in% snps$snp]
    tmp2=tmp[snp %in% snps$asnp]
    tmp2[,beta:=-beta][,snp:=paste(chromosome,pos38,ALT,REF,sep="-")]
    tmp=rbind(tmp1,tmp2)
    tmp[!duplicated(snp)]
})

variants=lapply(DATA, "[[", "snp")  %>% unlist()  %>% unique()

## the Chiou and UKBB datasets are too big. filter rare snps, not found in other datasets, that have large p > 0.01
ovariants=lapply(DATA[-c(6,10)], "[[", "snp")  %>% unlist()  %>% unique()
maf6=pmin(DATA[[6]]$effect_allele_frequency, 1 - DATA[[6]]$effect_allele_frequency)
maf10=pmin(DATA[[10]]$effect_allele_frequency, 1 - DATA[[10]]$effect_allele_frequency)
DATA[[6]][,P:=as.numeric(P)]
with(DATA[[6]], table(snp %in% ovariants | maf6 > 0.01 | P <0.01))
with(DATA[[10]], table(snp %in% ovariants | maf10 > 0.01 | P <0.01))
DATA[[6]]=DATA[[6]][snp %in% ovariants | maf6 > 0.01 | P <0.01]
DATA[[10]]=DATA[[10]][snp %in% ovariants | maf10 > 0.01 | P <0.01]
variants=lapply(DATA, "[[", "snp")  %>% unlist()  %>% unique()
length(variants)
save(DATA, file="gwas.RData")
