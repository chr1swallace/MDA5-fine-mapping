library(data.table)
library(magrittr)
## NB IFIH1 is at chr2:162,267,074-162,318,684 on build 38
st=161800000 # looking at manhattans, don't need a full 1 mb on left. this will make LD easier to calculate
en=163318684
message("window: ",st,"-",en)
source("get-data.R")

## * flip alleles
    ## Yorgo: I want to list the MDA5 alleles at residues 843 and 946 as T946A and R843H but in your meta analysis you used A946T and H843R. Can I just change the sign of beta to reverse the allele or does beta need to be recalculated?
v=fread("variants.txt")
   MR=lapply(DATA, function(x) 
        merge(x[,.(pos38,REF,ALT,beta,standard_error,P,disease,study)],
              v[,.(Variant.ID,SNP,pos38,ref_allele,alt_alleles)],
              by="pos38"))  %>%
        rbindlist()
    split(MR[,.(SNP,REF,ALT,study,beta,P)], MR$SNP)

MR[,snp:=paste(2,pos38,REF,ALT,sep="-")]
MR[,rsnp:=paste(2,pos38,ALT,REF,sep="-")]


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

toflip=c("2-162272314-T-C")
flipped=c("2-162272314-C-T")
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
## add splice donor sites and one tag
uMR=rbind(uMR,
          data.table(snp=c("2-162268086-C-T","2-162279995-C-G","2-162280432-A-G"),
                     rsnp=c("",""),
                     SNP=c("rs35732034", "rs35337543","rs72871627"),
                     pos38=c(162268086,162279995,162280432)))
uMR[,MAF:=ifelse(is.na(MAF[snp]), MAF[rsnp], MAF[snp])]
uMR[,EAF:=MAF]
uMR[,MAF:=pmin(MAF,1-MAF)] # don't do this, but call it EAF
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

## write out data we are about to plot
nm=names(D)
nm[names(D)=="liu_uc"]="Ulcerative colitis"
nm[names(D)=="liu_cro"]="Crohn's disease"
nm[names(D)=="robs_t1d_eur"]="Type 1 diabetes (Eur)"
nm[names(D)=="chiou_t1d"]="Type 1 diabetes"
nm[names(D)=="tsoi_pso"]="Psoriasis"
nm[names(D)=="ukbb_hypo"]="Hypothyroidism"

Duse=grep("liu|robs_t1d_eur|tsoi|hypo",names(D))
for(i in Duse) {
    dt=as.data.table(D[[i]][c("snp","position","beta","varbeta")])
    dt[,se:=sqrt(varbeta)][,varbeta:=NULL][,study:=nm[i]]
    fwrite(dt, file=paste0(names(D)[i],"-v2.csv"))
}

## * supplementary figure - manhattans

pdf("manhattans5.pdf",height=8,width=12,bg="white")
par(mfrow=c(2,3)) 
for(i in Duse) {
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

B=lapply(Duse, function(i) {
    d=D[[i]]
    keep=which(d$position %in% uMR$pos38)
    data.table(disease=names(D)[i], #nm[[i]],
               snp=d$snp[keep],
               position=d$position[keep],
               beta=d$beta[keep],
               se=sqrt(d$varbeta[keep]))
    })  %>% rbindlist()
B=merge(B, uMR[,.(position=pos38,SNP)])
dropsnps=c("rs35337543","R843H","H843R") # won't include in the Bonferroni, because these are tags of opther variants we will include
w=dcast(B[!(SNP %in% dropsnps & disease %in% names(D)[Duse])], position+ snp ~ disease, value.var=c("beta","se"))
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
## want to select one dataset per disease per snp. preferably the fine mapping dataset. otherwise the largest - so this is chiou_t1d when robs_t1d_eur is missing
levels=c(names(D)[Duse],"chiou_t1d","forgetta_t1d")
levels=c(levels,setdiff(names(D),levels))
B[,disease:=factor(disease,levels=levels)]
B=B[order(SNP,disease)]
B[,disease_short:=tstrsplit(disease,"_")[[2]]]
B[,totest:=!duplicated(disease_short),by="SNP"]
B[SNP %in% dropsnps,totest:=FALSE]
B[,.(SNP,disease,disease_short,totest)]
for(i in 1:nrow(B)) {
    B$inset[i]=B$snp[i] %in% unlist(lapply(S[[as.character(B$disease)[i]]]$sets$cs,names))
}
B[,EAF:=MAF[snp]]

## only Robertson EUR from t1d studies
B[,totest:=!(SNP %in% c("R843H","rs35337543"))]
B[totest==TRUE,padj:=pmin(1,p*nrow(B[totest==TRUE]))]
fwrite(B[totest==TRUE | SNP=="rs35337543" | (SNP=="R843H" & disease %in% B[totest==TRUE]$disease)][order(SNP,disease),.(SNP,snp,EAF,disease=disease_short,study=disease,totest,beta,se,inset,p,padj)], file="supptable-robertson.csv")



## * look at effect estimates
#B=B[SNP!="N160D"] # not significantly associated, but still interested
B[,dnm:=disease]
for(i in seq_along(D))
    B[disease==names(D)[i], dnm:=nm[i]]

m=melt(w, c("position","snp","beta_ukbb_hypo","se_ukbb_hypo"),
       list(beta=c("beta_liu_cro", "beta_liu_uc", "beta_tsoi_pso", "beta_robs_t1d_eur"),
            se=c("se_liu_cro", "se_liu_uc", "se_tsoi_pso", "se_robs_t1d_eur")))
m[,disease:=c("Crohn's disease","Ulcerative colitis", "Psoriasis", "Type 1 diabetes")[variable]]
#m=m[snp!="2-162310909-T-C"] # N160D is not significant anywhere
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

## what is linear coefficient?
library(pheatmap)
variants=unique(B[padj<0.05 & !(SNP %in% dropsnps)]$snp) # manually checked, already ordered by position
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

