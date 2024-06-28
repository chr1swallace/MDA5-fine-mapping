library(data.table)
library(magrittr)
## NB IFIH1 is at chr2:162,267,074-162,318,684 on build 38
st=161800000 # looking at manhattans, don't need a full 1 mb on left. this will make LD easier to calculate
en=163318684
message("window: ",st,"-",en)
source("get-data.R")

## * flip alleles
v=fread("variants.txt")
   MR=lapply(DATA, function(x) 
        merge(x[,.(pos38,REF,ALT,beta,standard_error,P,disease,study)],
              v[,.(Variant.ID,SNP,pos38,ref_allele,alt_alleles)],
              by="pos38"))  %>%
        rbindlist()
    split(MR[,.(SNP,REF,ALT,study,beta,P)], MR$SNP)

MR[,snp:=paste(2,pos38,REF,ALT,sep="-")]
MR[,rsnp:=paste(2,pos38,ALT,REF,sep="-")]

    ## Yorgo: I want to list the MDA5 alleles at residues 843 and 946 as T946A and R843H but in your meta analysis you used A946T and H843R. Can I just change the sign of beta to reverse the allele or does beta need to be recalculated?

    MR[SNP=="A946T",REF:="T"]
    MR[SNP=="A946T",ALT:="C"]
    MR[SNP=="A946T",beta:=-beta]
    MR[SNP=="A946T",SNP:="T946A"]
    MR[SNP=="H843R",REF:="C"]
    MR[SNP=="H843R",ALT:="T"]
    MR[SNP=="H843R",beta:=-beta]
    MR[SNP=="H843R",SNP:="R843H"]

toflip=c("2-162267541-C-T")
flipped=c("2-162267541-T-C")
DATA=lapply(DATA, function(d) {
    d[snp==toflip,c("beta","REF","ALT"):=list(-beta,ALT,REF)]
    d[snp==toflip,snp:=flipped]
    d
})
LD[,toflip]=-LD[,toflip]
LD[toflip,]=-LD[toflip,]
colnames(LD)  %<>% sub(toflip,flipped,.)
rownames(LD)  %<>% sub(toflip,flipped,.)
names(MAF)  %<>% sub(toflip,flipped,.)

uMR=MR[!duplicated(snp), .(snp,rsnp,SNP,pos38)]
## add splice donor sites
uMR=rbind(uMR,
          data.table(snp=c("2-162268086-C-T","2-162279995-C-G","2-162280432-A-G"),
                     rsnp=c("",""),
                     SNP=c("rs35732034", "rs35337543","rs72871627"),
                     pos38=c(162268086,162279995,162280432)))
uMR[,MAF:=ifelse(is.na(MAF[snp]), MAF[rsnp], MAF[snp])]
uMR[,MAF:=pmin(MAF,1-MAF)]
uMR[,description:=ifelse(grepl("rs",SNP),"splice donor","loss of function")]    
uMR[SNP=="rs72871627",description:="r2=0.99 with rs35337543"]
uMR

## make datasets
library(coloc)
make_dataset=function(d) {
    d=d[ snp %in% names(MAF) ]
    list(snp=d$snp,
         position=d$pos38,
         beta=d$beta,
         varbeta=d$standard_error^2,
         N=59e+3,s=0.4, # approx
         MAF=MAF[d$snp],
         LD=LD[d$snp,d$snp],
         type="cc")
}
D=lapply(names(DATA), function(nm) { message(nm); make_dataset(DATA[[nm]])})
## names(D)=grep("t1d|uc|cro",names(DATA),value=TRUE)
names(D)=names(DATA)
save(D,file="~/scratch/ifih1/D.RData")

s=runsusie(D[[7]])
s$sets$cs

S=lapply(D, runsusie)
save(S,file="~/scratch/ifih1/S.RData")

save(LD,MAF,MR,uMR,variants,file="extra.RData")

## * examine results
(load("~/scratch/ifih1/S.RData"))
(load("~/scratch/ifih1/D.RData"))
(load("extra.RData")) # flipped LD, MAF in here

## compare credsets
sapply(S, function(s) length(s$sets$cs))
sapply(S, function(s) sapply(s$sets$cs,names))

CS=lapply(S, function(s) s$sets$cs)
POS=lapply(CS, function(cs) lapply(cs, function(x) as.numeric(gsub("^2-|-[AGCT]-[ACGT]","",names(x)))))
VAR=lapply(POS, function(pos) lapply(pos, function(p) v[pos38 %in% p]))

library(coloc)
# compare everything to t1d_eur as strongest signals and most credsets
a=coloc.susie(S$liu_uc,S[[7]]) #$summary
a$summary[PP.H4.abf>.8] # 1-1 2-2 3-3

a=coloc.susie(S$tsoi_pso,S[[7]]) #$summary
a$summary[PP.H4.abf>.8] # 1-1 

a=coloc.susie(S$ukbb_hypo,S[[7]]) #$summary
a$summary[PP.H4.abf>.8] # 1-1 2-2 3-3

## what is 3rd (shared) t1d credset and 4th (unique) t1d cs?
S[[7]]$sets$cs

## $L1
## 2-162267541-C-T  = T946A
##             591 

## $L4
## 2-162234023-G-T 2-162238769-C-T 2-162268086-C-T - splice donor variant (3rd snp)
##             463             485             595 

## $L2
## 2-162279995-C-G 2-162280432-A-G - splice donor variant (1st snp)
##             639             643 

## $L3
## 2-162268127-T-C 2-162275836-C-T  - I923V (1st snp)
##             597             627


## t1d
a=coloc.susie(S[[5]],S[[8]]) # S5 has one set, matches set 1 for S8
a=coloc.susie(S[[6]],S[[8]]) # S6 has 3 sets. matches S6-S8 are 1-1 3-2 2-3
a=coloc.susie(S[[7]],S[[8]]) # S7 has 4 sets. matches S7-S8 are 1-1 2-2 3-3 4-4

## t1d-uc
a=coloc.susie(S[[1]],S[[8]]) # S1 has 3 sets. matches S1-S8 are 1-1 2-2 3-3
a=coloc.susie(S[[3]],S[[8]]) # S3 has 2 sets. no matches with S8
a$summary[PP.H4.abf>.8] 

## t1d-pso
a=coloc.susie(S[[8]],S[[9]]) # pso has 1 set. matches 1-1
a$summary[PP.H4.abf>.8] 

## t1d-hypo
a=coloc.susie(S[[8]],S[[10]]) # hypo has 3 sets. matches to S8 are 1-1 2-2 3-3
a$summary[PP.H4.abf>.8] 

cols=c("dodgerblue2", 
        "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2", 
        "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", #"gray70", 
        "khaki2", "maroon", "orchid1", "deeppink1", "blue1", 
        "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", 
       "darkorange4", "brown")
cols[8]=cols[7] # because of the tag of splice variant rs35337543

library(grid)
library(gridBase)
library(gridExtra)

## look at key snps from Yorgo
v=fread("variants.txt")
 
pdf("manhattans5.pdf",height=8,width=12,bg="white")
nm=names(D)
nm[names(D)=="liu_uc"]="Ulcerative colitis"
nm[names(D)=="liu_cro"]="Crohn's disease"
nm[names(D)=="robs_t1d_eur"]="Type 1 diabetes"
nm[names(D)=="tsoi_pso"]="Psoriasis"
nm[names(D)=="ukbb_hypo"]="Hypothyroidism"
par(mfrow=c(2,3)) 
for(i in grep("liu|robs_t1d_eur|tsoi|hypo",names(D))) {
    plot_dataset(D[[i]],main=nm[i])
    p=2*pnorm(-abs(D[[i]]$beta)/sqrt(D[[i]]$varbeta))
    pos=D[[i]]$position
    ## credible sets
    highlight_list = lapply(S[[i]]$sets$cs, names)
    for(j in seq_along(highlight_list)) {
        l=highlight_list[[j]]
        wl=which(uMR$snp %in% l | uMR$rsnp %in% l)
        if(length(wl)>1) { wl=wl[1] }
        colj=if(length(wl)) { cols[wl] } else { cols[10 + j] }
        mm=match(l, D[[i]]$snp) # all snps in the set
        points(pos[mm], -log10(p[mm]), col=colj, pch="o", cex=2) 
        if(length(wl)) {
            mm=match(uMR$pos38[wl], pos) # missense snp in the set
            text(pos[mm],-log10(p[mm]),label=uMR$SNP[wl],pos=4,col=colj)
        }
    }
    ## missense
    w=which(pos %in% uMR$pos38[1:7])
    points(pos[w],-log10(p[w]),pch="x")
    mm=match(pos[w],uMR$pos38)
    if(!length(highlight_list)) {
        text(pos[w],-log10(p[w]),label=uMR$SNP[mm],pos=4)
        next
    }
    for(j in seq_along(w)) {
        snpj=D[[i]]$snp[w[j]]
        inset=sapply(highlight_list, function(l) snpj %in% l)  %>% any()
        if(!inset)
            text(D[[i]]$position[w[j]],-log10(p[w[j]]),label=uMR$SNP[mm[j]],pos=4)
    }
}

frame()
# Grid regions of current base plot (ie from frame)
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
# Table grob
grob <-  uMR[SNP %in% c("T946A","I923V","rs35337543","rs72871627","rs35732034"),
             .(SNP,MAF,description)]  %>% as.data.frame()  %>% tableGrob(.,rows=NULL)
grid.draw(grob)

dev.off()

## ** plot region

library(Gviz)
library(biomaRt)
bm <- useEnsembl(#host = "https://grch38.ensembl.org", 
              biomart = "ENSEMBL_MART_ENSEMBL", 
              dataset = "hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome = "hg38",
                                    #chromosome = 'chr2', start = st, end = en,
                                    symbol="IFIH1", name = "ENSEMBL", biomart = bm, stacking="squish",
                                    filter=list(hgnc_symbol="IFIH1"))
aTrack <- AnnotationTrack(start = uMR$pos38[1:7], width = 1, 
                          chromosome = "chr2", 
                          strand = c("*","*","*","*","*","*","*"),
                          id = uMR$SNP[1:7],
                          genome = "hg38", name = "SNPs",
                          group=uMR$SNP[1:7])
feature(aTrack)=group(aTrack)
## plotTracks(aTrack)
## plotTrcks(biomTrack)
pdf("gene_snps.pdf",height=4,width=12/3,bg="white")
plotTracks(list(biomTrack, aTrack), transcriptAnnotation = "symbol", groupAnnotation="group",
           T964A=cols[1], I923V=cols[2],rs35732034=cols[6], rs35337543=cols[7] )
dev.off()

## ** Bonferroni analysis of missense + splice donors 

B=lapply(c(1,2,7,9,10), function(i) {
    d=D[[i]]
    keep=which(d$position %in% uMR$pos38)
    data.table(disease=names(D)[i], #nm[[i]],
               snp=d$snp[keep],
               position=d$position[keep],
               beta=d$beta[keep],
               se=sqrt(d$varbeta[keep]))
    })  %>% rbindlist()
B=merge(B, uMR[,.(position=pos38,SNP)])
B=B[SNP!="rs35337543"]
B=B[SNP!="R843H"]
B=B[SNP!="H843R"]
w=dcast(B, position+ snp ~ disease, value.var=c("beta","se"))
m=melt(w, c("position","snp","beta_robs_t1d_eur","se_robs_t1d_eur"),
       list("beta_liu_cro", "se_liu_cro",
            "beta_liu_uc", "se_liu_uc",
            "beta_tsoi_pso", "se_tsoi_pso",
            "beta_ukbb_hypo",       "se_ukbb_hypo"))
setnames(m, paste0("value",1:8), c("beta_liu_cro", "se_liu_cro",
            "beta_liu_uc", "se_liu_uc",
            "beta_tsoi_pso", "se_tsoi_pso",
            "beta_ukbb_hypo",       "se_ukbb_hypo"))
m=merge(m, uMR[,.(position=pos38,SNP)])
## fix effect for rs35337543 from co-credset snp
#toput=m[SNP=="rs35337543 *", 6:11, with=FALSE]
B[,z:=beta/se][,p:=2*pnorm(-abs(z))]
B[,padj:=pmin(1,p*nrow(B))]
for(i in which(B$disease!="liu_cro"))
    B$inset[i]=B$snp[i] %in% unlist(lapply(S[[B$disease[i]]]$sets$cs,names))
options(width=100)
B[order(padj),.(SNP,snp,disease,beta,inset,p,padj)]

## try again
B=B[SNP!="N160D"]
B[,dnm:=disease]
for(i in seq_along(D))
    B[disease==names(D)[i], dnm:=nm[i]]

m=melt(w, c("position","snp","beta_ukbb_hypo","se_ukbb_hypo"),
       list(beta=c("beta_liu_cro", "beta_liu_uc", "beta_tsoi_pso", "beta_robs_t1d_eur"),
            se=c("se_liu_cro", "se_liu_uc", "se_tsoi_pso", "se_robs_t1d_eur")))
m[,disease:=c("Crohn's disease","Ulcerative colitis", "Psoriasis", "Type 1 diabetes")[variable]]
m=m[snp!="2-162310909-T-C"] # N160D is not significant anywhere
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())

## pick a common scale
dcols=cols[1:5]
names(dcols)=unique(B$dnm)

plin=ggplot(m, aes(x=beta_ukbb_hypo,y=beta, col=disease)) + geom_point() +
    geom_linerange(aes(ymin=beta-1.96*se,ymax=beta+1.96*se)) +
    ## geom_linerange(aes(xmin=beta_robs_t1d_eur-1.96*se_robs_t1d_eur,xmax=beta_robs_t1d_eur+1.96*se_robs_t1d_eur)) +
    background_grid() +
    geom_abline(linetype="dashed") +
    geom_abline(slope=-1,linetype="dashed") +
    geom_smooth(method="lm",se=FALSE) +
    scale_colour_manual(values=dcols, breaks=names(dcols)) +
    labs(x="log OR (HypoThy)",y="log OR (disease)")
    ## ylim(-.7,0.5)



## TODO estimate slope for comparing diseases using MR IVW

## what is LD between these snps?
library(pheatmap)
variants=unique(B[padj<0.05]$snp) # manually checked, already ordered by position
diseases=unique(B$disease)
M=matrix(1,length(diseases),length(diseases))
dimnames(M)=list(diseases,diseases)
P=M
todo=expand.grid(d1=diseases,d2=diseases) %>% as.data.table()
todo=todo[d1!=d2]

library(colocPropTest)
for(i in 1:nrow(todo)) {
    d1=todo$d1[i]
    d2=todo$d2[i]
    m=merge(B[disease==d1,.(SNP,snp,beta,se)],
            B[disease==d2,.(SNP,snp,beta,se)], by=c("SNP","snp"), suffixes=c(1,2))
    V1=m$se1 %*% t(m$se1) * LD[m$snp,m$snp]
    V2=m$se2 %*% t(m$se2) * LD[m$snp,m$snp]
    res=colocPropTest:::estprop_slow(m$beta1, m$beta2, V1, V2)$result
    M[d1,d2]=res["eta.hat"]
    P[d1,d2]=pchisq(res["chisquare"], df=res["n"]-1)
}

summary(as.vector(P)) # all are proportional
rownames(M)  %<>% sub("liu_uc","UC",.)
rownames(M)  %<>% sub("liu_cro","CRO",.)
rownames(M)  %<>% sub("robs_t1d_eur","T1D",.)
rownames(M)  %<>% sub("tsoi_pso","Pso",.)
rownames(M)  %<>% sub("ukbb_hypo","HThy",.)
colnames(M)  %<>% sub("liu_uc","UC",.)
colnames(M)  %<>% sub("liu_cro","CRO",.)
colnames(M)  %<>% sub("robs_t1d_eur","T1D",.)
colnames(M)  %<>% sub("tsoi_pso","Pso",.)
colnames(M)  %<>% sub("ukbb_hypo","HThy",.)
## for(i in seq_along(D)) {
##     colnames(M) %<>% sub(names(D)[i], nm[i], .)
##     rownames(M) %<>% sub(names(D)[i], nm[i], .)
## }

library(RColorBrewer)
paletteLength <- 50
myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- seq(-4.5,4.5,length.out=50)
grob=pheatmap(M,display_numbers=TRUE,fontsize=12,number_color="grey20",fontsize_number=12,
         cluster_rows=FALSE,cluster_cols=FALSE,
         color=myColor, breaks=myBreaks, legend=FALSE)

## now plot together with effect sizes
pdf("effects-by-snp.pdf",bg="white",height=8,width=10)
ggplot(B[SNP %in% B[padj<0.05]$SNP], aes(y=dnm,x=beta,xmin=beta-1.96*se,xmax=beta+1.96*se,col=dnm)) + geom_pointrange() + facet_wrap(~SNP) + #,scales="free_x") +
    geom_vline(xintercept=0) +
    scale_colour_manual("Disease",values=dcols,breaks=names(dcols),guide = guide_legend(reverse = TRUE)) +
    labs(x="log odds ratio",y="Disease") +
    background_grid(major="y") +
    theme(legend.position="none") 
          ## legend.justification="bottom",
          ## legend.background = element_blank(),
          ## legend.box.background = element_rect(colour = "black"))
# Grid regions of current base plot (ie from frame)
vp <- viewport(x=2/3+.05 ,y=0.1,width=1/3-.05, height=1/2-.2,just=c("left","bottom"))
pushViewport(vp)
grid.draw(as_grob(plin + theme(legend.position="none")))
dev.off()

## * conditional fine mapping with gcta
## R843H and I923V appear pretty high on manhattans if typed, but not in credible sets. is this just due to LD? use cojo to check.

write_cojo=function(d) {
##    ID      ALLELE1 ALLELE0 A1FREQ  BETA    SE      P       N
## chr1:11171:CCTTG:C      C       CCTTG   0.0831407       -0.0459889      0.0710074       0.5172  180590
## chr1:13024:G:A  A       G       1.63957e-05     -3.2714 3.26302 0.3161  180590 
    ss=tstrsplit(d$snp,"-")
    z=d$beta/sqrt(d$varbeta)
    dt=with(d,data.table(ID=snp,ALLELE1=ss[[3]],ALLELE2=ss[[4]],A1FREQ=MAF[ snp ], BETA=beta, SE=sqrt(varbeta), P=pnorm(-abs(z))*2, N=N))
    tmp=tempfile()
    fwrite(dt, file=tmp, sep="\t")
    return(tmp)
}

## fix up plink files
snps[,zero:=0]
fwrite(snps[,.(chr,snp,zero,position,a0,a1)],
       file="ukbb40k_alt.bim",
       sep="\t", col.names=FALSE)
if(!file.exists("ukbb40k_alt.bed"))
    paste0("ln -s ukbb40k.bed ukbb40k_alt.bed")  %>% system()
if(!file.exists("ukbb40k_alt.fam"))
    paste0("ln -s ukbb40k.fam ukbb40k_alt.fam")  %>% system()

    
stub="ukbb40k_alt"
for(nm in names(D)) {
    message("\nRunning cojo for ",nm)
    d=D[[nm]]
    tmp=write_cojo(D[[nm]])
    tmpout=tempfile()
    paste0("gcta-1.94.1 --cojo-file ",tmp, #" --chr ",chr,
           " --bfile ",stub,
           " --maf 0.005 --cojo-slct --cojo-p 1e-3 --out ",tmpout)  %>% system()
    cojo_tmpout=paste0(tmpout,".jma.cojo")
    if(!file.exists(cojo_tmpout)) {
        cleanup(tmpout,tmp)
        next
    }
    csnps=fread(cojo_tmpout)

    if(nrow(csnps)==1) {
        COND[[nm]]=list(with(d,data.table(snp=snp,beta=beta,varbeta=varbeta)))
        attr(COND[[nm]],"condsnps")=NULL
    } else {
        tmpsnps=tempfile()
        ## HERE - conditional snp will be missing - need to deal with that in format of output
        cresults=lapply(1:nrow(csnps), function(i) { 
            fwrite(csnps[-i,.(SNP)],file=tmpsnps,col.names=FALSE)
            paste0("gcta-1.94.1 --cojo-file ",tmp, #" --chr ",chr,
                   " --bfile data/reference/byblock/",args$block,"_ukbb40k_alt_hg38",
                   " --maf 0.005 --cojo-cond ",tmpsnps," --out ",tmpout)  %>% system()
            result=paste0(tmpout,".cma.cojo")  %>% fread()
            result[,.(snp=SNP,beta=bC,varbeta=bC_se^2)]
        })  
        attr(cresults,"condsnps")=csnps
        COND[[nm]]=cresults
        unlink(tmpsnps)
    }
    cleanup(tmpout,tmp)
}
            
save(COND,file=outfile_cojo)

message("COJO complete")



## * comparison between traits
## uc
library(coloc)
lapply(S[c(1,3)], function(s) s$sets$cs) # no overlap. check with coloc
coloc.abf(D[[1]],D[[3]]) ## pp.h4=0.959
a=coloc.susie(S[[1]],S[[3]]) #$summary
a$summary # liu signal 1 matches delange signal 1. no other matches.

## t1d
a=coloc.susie(S[[5]],S[[8]]) # S5 has one set, matches set 1 for S8
a=coloc.susie(S[[6]],S[[8]]) # S6 has 3 sets. matches S6-S8 are 1-1 3-2 2-3
a=coloc.susie(S[[7]],S[[8]]) # S7 has 4 sets. matches S7-S8 are 1-1 2-2 3-3 4-4

## t1d-uc
a=coloc.susie(S[[1]],S[[8]]) # S1 has 3 sets. matches S1-S8 are 1-1 2-2 3-3
a=coloc.susie(S[[3]],S[[8]]) # S3 has 2 sets. no matches with S8
a$summary[PP.H4.abf>.8] 

## t1d-pso
a=coloc.susie(S[[8]],S[[9]]) # pso has 1 set. matches 1-1
a$summary[PP.H4.abf>.8] 

## t1d-hypo
a=coloc.susie(S[[8]],S[[10]]) # hypo has 3 sets. matches to S8 are 1-1 2-2 3-3
a$summary[PP.H4.abf>.8] 

plot_dataset(D[[8]],S_orig[[1]]) # v1
plot_dataset(D[[8]],S_orig[[2]]) # v1 eur
plot_dataset(D[[8]],S[[8]]) # v2
plot_dataset(D[[8]],S[[7]]) # v2 eurS


```

## Results

Plot data + susie results.  Credible sets are coloured and numbered.  SNPs nominated by Rahul are plotted as black points and labelled above by "R".  Lead SNPs in the Robertson paper are plotted as "X".
It looks like the fine mapping failed for Chiou (mismatch in LD or data quality?) and is weak in Forgetta (presumably sample size). 

In the Robertson results, all 4 lead SNPs from Robertson are captured + 1 extra, which is a Rahul nomination.  The results in the full Robertson data and the European subset are in broad agreement, which is comforting (the LD we use is from a European reference panel, so important to check that this doesn't introduce errors).

```{r,include=TRUE}
cs_snps=lapply(S, function(s) lapply(s$sets$cs,names))

library(coloc)
par(mfrow=c(2,2))
for(i in 1:4) {
plot_dataset(D[[i]],S[[i]])
m=match(lv$snp38,D[[i]]$snp) %>% setdiff(.,NA)
p=with(D[[i]], pnorm(abs(beta/sqrt(varbeta)),lower=FALSE)*2)
text(x=D[[i]]$position[m],y=-log10(p[m]),label="R",pos=3) # 4/4 overlap - remarkable
points(x=D[[i]]$position[m],y=-log10(p[m]),pch=16) # 4/4 overlap - remarkable
m=match(rob.leads$snp38,D[[i]]$snp) %>% setdiff(.,NA)
points(x=D[[i]]$position[m],y=-log10(p[m]),pch="X") # 4/4 overlap - remarkable
title(main=names(D)[i])
}
```

Check PP

``` r
for(i in 1:4) {
    message(names(S)[i])
    print(S[[i]]$pip[ unlist(lapply(S[[i]]$sets$cs, names)) ])
}
```

```{r,include=TRUE,fig.show="asis"}
for(i in 1:4) {
  inn=paste0("in.",names(D)[i])
  csn=paste0("cs.",names(D)[i])
  pn=paste0("p.",names(D)[i])
  bn=paste0("beta.",names(D)[i])
  b=D[[i]]$beta
  z=b/sqrt(D[[i]]$varbeta)
  p=pnorm(-abs(z))*2
  names(b)=names(p)=D[[i]]$snp
  lv[[inn]]=lv$snp38 %in% D[[i]]$snp
  lv[[csn]]=lv$snp38 %in% unlist(cs_snps[[i]])
  lv[[bn]]=b[lv$snp38]
  lv[[pn]]=p[lv$snp38]
}

maf=rowMeans(haps,na.rm=TRUE)
lv$MAF.1kg=maf[lv$snp38]
results=lv[,c("SNP","Variant.ID","snp38","MAF.1kg",#"MDA5 Effect",
              grep("^in.|^cs.|^p.|^b.",names(lv),value=TRUE)),with=FALSE]
results=results[, !grepl("_eur",names(results)),with=FALSE]
kable(results)
fwrite(results, file="ifih1-t1d-results.csv")

## with(rb[snp %in% v$snp38])

```

Check LD between key snps and Robertson cs snps

``` r
thaps=t(haps[unique(intersect(c(lv$snp38,unlist(cs_snps[[3]])),rownames(haps))),])
sdsnps=apply(thaps,2,sd,na.rm=TRUE)
if(any(sdsnps==0)) {
  w=which(sdsnps==0 | is.na(sdsnps))
  message("monomorphic: ",colnames(thaps)[w])
  thaps=thaps[,-w]
}
ld.rob=cor(thaps,use="pair")
ld.rob[,"2-162310909-T-C"]^2 # all < 0.01

```

Examine Robertson results at the fine mapped signals. "pip" is the posterior probability of inclusion for each variant in the set. Higher=greater likelihood of being causal.  Of the SNPs in set 4, which don't contain any nominated by Rahul, only 2:162268086 (rs35732034) is in IFIH1, and is a splice donor according to dbSNP [[https://www.ncbi.nlm.nih.gov/snp/rs35732034]].  

```{r,include=TRUE,results="asis"}
sets=lapply(S$Robertson$sets$cs,names)
it=1
for(idx in seq_along(S$Robertson$sets)) {
  snps=names(S$Robertson$sets$cs[[idx]])
  pos=as.numeric(tstrsplit(snps,"-")[[2]])
  tmppip=S$Robertson$alpha[ S$Robertson$sets$cs_index[idx], snps]
  srobs=robs[pos38 %in% pos]
  srobs=merge( v[,.(pos38,protein)], srobs[,.(pos38,effect_allele,other_allele, beta, p_value)],by="pos38",all.y=TRUE)[match(pos,pos38)]
  srobs[,pip:=tmppip]
  print(kable(srobs,
              caption=paste0("Cred set ",it)))
  cat("\n")
  it=it+1
}
```

## Imputation of missing variants

Consider summary statistic imputation of 2-162277580-C-A which is missing from Robertson.

There are only 5 heterozygotes in our reference panel

```{r,include=TRUE}
table(haps[v$snp38[3],])
```

But there are several SNPs in LD with it
```{r,include=TRUE}
thaps=t(haps)
hsd=apply(thaps,2,sd,na.rm=TRUE)
thaps=thaps[,hsd>0]

tld=cor(thaps[,-which(colnames(thaps)==v$snp38[3])],thaps[,v$snp38[3]])[,1]
print(tld[tld>0.5])
```

whilst some of these are in Robertson
```{r,include=TRUE}
d=S[["Robertson"]]
print(names(tld)[tld>0.5] %in% d$snp)
```
none make it into any credible sets
```{r,include=TRUE}
s=S$Robertson
cs_snps=lapply(s$sets$cs,names) 
lapply(cs_snps, intersect, names(tld)[tld>0.5])
```

Solution: use larger reference panel.  I have applied to access the Haplotype Reference Consortium data.

<!-- ``` -->
<!-- ## TODO: -->
<!-- ## 1. other diseases -->
<!-- ## lupus 12k EAS -->
<!-- if(!file.exists("33536424-GCST90011866-EFO_0002690.h.tsv.gz")) -->
<!--   system("wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90011001-GCST90012000/GCST90011866/harmonised/33536424-GCST90011866-EFO_0002690.h.tsv.gz") -->
<!-- lupus.wang=fread("33536424-GCST90011866-EFO_0002690.h.tsv.gz") -->
<!-- lupus.wang[,disease:="SLE"][,study:="Wang"] -->
<!-- M$lupus.wang=merge(v,lupus.wang[hm_chrom==2,.(hm_pos,hm_variant_id,beta=hm_beta,P=p_value,REF=other_allele,ALT=effect_allele,disease,study)],by.y="hm_pos",by.x="pos38") -->
<!-- ## langefield lupus immunochip 18k -->
<!-- lupus.langefield=fread("~/share/Data/GWAS-summary/sle-immunochip/ea.imputed.out.gz") -->
<!-- ## bentham EUR gwas -->
<!-- lupus.bentham=fread("~/share2/02-Processed/SLE_Bentham_26502338_1-hg38.tsv.gz") -->
<!-- lupus.bentham[,disease:="SLE"][,study:="Bentham"] -->
<!-- M$lupus.bentham=merge(v,lupus.bentham[hm_CHR==2,.(hm_BP,hm_variant_id,beta=hm_BETA,P=P,REF,ALT,disease,study)],by.y="hm_BP",by.x="pos38") -->

<!-- Mr=rbindlist(M,fill=TRUE) -->
<!-- options(width=120) -->
<!-- Mr[order(protein),.(protein,snp38,REF,ALT,disease,study,beta,P)] -->

<!-- ``` -->

## Notes from meeting 12/8/2021

there are a pair of variants thought to work together. Rahul will send paper, I will check in ImmunoChip data

there are 3 more variants from Rahul to come

I will wait for HRC data and try ssimp

Rahul will check directionality of DNA risk on T1D vs predicted function at protein level.

## COVID variants

build 38
```{r}
library(data.table)
X=lapply(c("A2","B1","B2"), function(nm) {
    f=paste0("~/share/Data/GWAS-summary/covid19hg/COVID19_HGI_",nm,"_ALL_leave_23andme_20220403.tsv.gz")
    command=paste0("zcat ",f," | awk '$1==2 && $2 > ", min(v$pos38)-1000," && $2 < ",max(v$pos38) + 1000, "'")
    x=fread(cmd = command)
    x[,snp38:=gsub(":","-",V5)]
    x[ V2 %in% lv$pos38 | V18 %in% lv$Variant.ID ]
})
names(X)=c("A2","B1","B2")
lapply(X, function(x) merge(x[,.(snp38,b=V7,p=V9,maf=V17,rsid=V18)],
                            lv[,.(snp38,SNP)],by="snp38"))



merge(x[ V2 %in% v$pos38 , .(pos38=V2, snp38,b=V7,p=V9,maf=V17,rsid=V18)], lv[,snp38,SNP], by="snp38")
merge(x[ V2 %in% lv$pos38 , .(pos38=V2, snp38,b=V7,p=V9,maf=V17,rsid=V18)], lv[,snp38,SNP], by="snp38")

```


From Brian: what about these variants in DHX58?
DHX58	chr17	42105987	G	A	missense		NA	.	DHX58:NM_024119.3:exon9:c.1000C>T:p.(Arg334Cys)	0.00459265	0.0228	.	.	0.0012	0.0001	0.014	0.00009487	rs76998797	.	22.8	3.313	high	0.951	possibly damaging	0.239	high	0.004	deleterious	0.243	high	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0|1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0|1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0|1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0|1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0|1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0|1	0	0|1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	yes
DHX58	chr17	42101837	C	T	missense		NA	.	DHX58:NM_024119.3:exon14:c.1961G>A:p.(Arg654His)	0.00259585	0.005	0.0029	0.002	0.004	0.0019	0.0016	0.0048	rs117253612	Score=312;Name=lod=20	23.1	3.313	high	0.005	benign	0.239	low	0.04	deleterious	0.243	high	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0|1	0	0	0	0	0	0	0	0	0	0	0	0|1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0|1	0	0	0	0	yes

```{r}
library(data.table)
f="~/share/Data/GWAS-summary/covid19hg/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz"
command=paste0("zcat ",f," | awk 'NR==1'")
nm=fread(cmd = command)
command=paste0("zcat ",f," | awk '$1==17 && $2 > 42101000 && $2 < 42106000'")
x=fread(cmd = command)
setnames(x, names(nm))
x[POS %in% c(42105987,42101837)]
```


Now check IBD

``` r
blocks=fread("~/E/blocks.txt")
b=blocks[V2==2 & V3 < 162267541 & V4 > 162267541 ]$V1# chr2_block93
f="~/E/data/byblock"
list.files(f,pattern=b)
(load(file.path(f,paste0("datasusie_",b,".RData"))))
names(SUSIE.FM)
w=names(SUSIE.FM)[c(10,26)]
SUSIE.FM[[w[1]]]$sets # IBD
SUSIE.FM[[w[2]]]$sets # UC
    ```
