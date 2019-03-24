### Scripts used for G4Hunter related data in Prorok et al, 2019
### by L. Lacroix, laurent.lacroix@inserm.fr
### 20190323

# You NEED TO CHANGE THE DEFAULT PATH OF THE SCRIPT ACCORDING TO YOUR SETUP
localpath <- 'puthereyourpath'

setwd(localpath)

# The whole script requires the following files in your working directory:
# 'mm10G4H2.RData' ; 'rd_ori.RData' ; 'mm10_QP37.RData'
# the SNS peak and bam files have been downloaded from http://rsat-tagc.univ-mrs.fr/g4/g4_data.html
# you also need to create a directory 'Script_Results/' in your working directory.
# Data will be saved in this directory

### Main functions #############################################################
source('./G4Hunter_mm10_Ori_RFunction.r')

#### Main Script ###############################################################
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
chnb=21
seqinf <- Seqinfo(seqnames=seqnames(genome)[1:chnb],seqlengths=seqlengths(genome)[1:chnb],isCircular=isCircular(genome)[1:chnb],genome='mm10')

### Import SNS_peaks data
# data were imported from http://rsat-tagc.univ-mrs.fr/g4/g4_data.html
RawReads <- read.table('./ori_macs2-sicer_DESEQ2_normReadsNumber_sortByFC.tab',header=T,sep='\t',as.is=T)
ori_reads<- rawread2gr(RawReads)
rm(RawReads)
# rename the enforced group in to reinforced (as in the paper)
ori_reads[ori_reads$ori_group=='enforced']$ori_group <- 'reinforced'
ori_reads$ori_group <- factor(ori_reads$ori_group, levels=c('insensitive','new','reinforced','reduced','suppressed'))

# Data cleaning (3 suppressed ori have a positive FC!)
ori_reads[ori_reads$ori_group=='suppressed' & ori_reads$ori_FC>0]
#GRanges object with 3 ranges and 4 metadata columns:
#      seqnames            ranges strand |  ori_group    ori_FC ori_DMSreads ori_PhDCreads
#         <Rle>         <IRanges>  <Rle> |   <factor> <numeric>    <numeric>     <numeric>
#  [1]    chr12 14659372-14660262      * | suppressed      4.03  2.900398762   43.68155229
#  [2]     chr6 69099246-69099530      * | suppressed      3.74  1.346158477   15.60864788
#  [3]     chr1 88256084-88256838      * | suppressed      1.76  166.5724552   469.9468123
# le'ts remove this
ori_reads <- ori_reads[!(ori_reads$ori_group=='suppressed' & ori_reads$ori_FC>0)]

ori_c <- resize(ori_reads,fix='center', width=1)
ori_list_c <- split(ori_c, ori_reads$ori_group)
Ori_rep <- round(elementNROWS(ori_list_c)/length(ori_c),3)*100
###

### G4Hunt
# beware if your setup is not multicore friendly
# G4H2 <- do.call(c,mclapply(1:chnb, function(i) mG4huntref(i,hl=2,with.seq=F),mc.cores=6))
# can be quite long depending on the number of CPU (20' for me)
# save(G4H2, file='mm10G4H2.RData')

# Possibility to load precomputed G4Hunter results
load('mm10G4H2.RData')

G4data <- G4H2
seqlevels(G4data, pruning.mode='coarse') <- seqlevels(seqinf)[1:chnb]
###

### Fig2B
ori1GR <- GRanges('chr11',ranges=IRanges(start=60139550,end=60140000),strand='+')
ori1d1.1 <- GRanges('chr11',ranges=IRanges(start=60139550,end=60139767),strand='+')
ori1d1.2 <- GRanges('chr11',ranges=IRanges(start=60139824,end=60140000),strand='+')
ori1d2.1 <- GRanges('chr11',ranges=IRanges(start=60139550,end=60139767),strand='+')
ori1d2.2 <- GRanges('chr11',ranges=IRanges(start=60139869,end=60140000),strand='+')

test.seq1 <- as.character(Views(genome,ori1GR))
test.seq2 <- paste0(as.character(Views(genome,ori1d1.1)),as.character(Views(genome,ori1d1.2)))
test.seq3 <- paste0(as.character(Views(genome,ori1d2.1)),as.character(Views(genome,ori1d2.2)))

seq1 <- reverseComplement(DNAString(test.seq1))
seq2 <- reverseComplement(DNAString(test.seq2))
seq3 <- reverseComplement(DNAString(test.seq3))
sco1 <- runmean(G4translate(seq1),k=25)
sco2 <- runmean(G4translate(seq2),k=25)
sco3 <- runmean(G4translate(seq3),k=25)
pdf('./Script_Results/Fig2B.pdf',width=10)
plot(sco1[130:280],type='l',ylim=c(-4,4),xlab='Position on chr11',ylab='G4Hscore',cex.axis=1.4,cex.lab=1.4,col='red',lwd=2,x=(60139870:60139720),xlim=c(60139870,60139720),font.lab=2,font=2,xaxs='i',yaxs='i')
lines(sco2[130:280],lty=2,col='blue',lwd=2,x=(60139870:60139720))
lines(sco3[130:280],lty=5,col='green4',lwd=2,x=(60139870:60139720))
abline(h=-2,lty=4)
abline(h=2,lty=4)
abline(h=0)
# Grey zone (-1:1) and (-1.5:1.5) in G4Hscore
polygon(c(60139870,60139720,60139720,60139870),c(1.5,1.5,-1.5,-1.5),col=rgb(0,0,0,0.1),border=F)
polygon(c(60139870,60139720,60139720,60139870),c(1,1,-1,-1),col=rgb(0,0,0,0.1),border=F)
# Deleted regions
lines(c(60139823,60139768),c(3.7,3.7),col=rgb(1,0,0,0.9),lwd=3,lend=2)
lines(c(60139823,60139768),c(3.7,3.7),col=rgb(0,0,1,0.8),lwd=3,lty=2,lend=2)
lines(c(60139868,60139768),c(3.4,3.4),col=rgb(1,0,0,0.9),lwd=3,lend=2)
lines(c(60139868,60139768),c(3.4,3.4),col=rgb(0,0.8,0.2,0.8),lwd=3,lty=5,lend=2)
legend('topright', legend=c('Ori1','Ori1 allele1','Ori1 allele2'), col=c('red','blue','green4'),bty='n',lty=c(1,2,5),lwd=2,cex=1.4)
dev.off()

###

### Fig3E
pdf('./Script_Results/Fig3E.pdf',width=21)
ylab1 <- c('Reads per origin','','','','')
par(mfcol=c(1,5))
bp <- lapply(c(1,2,3,4,5), function(i) {boxplot(ori_list_c[[i]]$ori_DMSreads,ori_list_c[[i]]$ori_PhDCreads,outline = F,ylim=c(0,900),main=paste0(toupper(names(ori_list_c)[i]),' (',Ori_rep[i],'%)'),axis=F,cex.axis=1.6,cex.main=2.5,font.axis=2);mtext(ylab1[i],2,2.5,font=2,cex=1.4);mtext(c('-','+ PhenDC3'),1,1,at=c(1,2),font=2,cex=1.2);axis(2,font.axis=2,cex.axis=1.6)})
dev.off()
###

### Fig4B
# First idenfiy Gaps of more than 20N in the genome to avoid shuffling if then
# beware if your setup is not multicore friendly

listGaps <- mclapply(1:chnb, function(i) {ra= findNgaps(genome[[i]]);resu=GRanges(seqnames=seqnames(genome)[i],ranges=ra,seqinfo=seqinfo(genome));return(resu)},mc.cores=4)
Ngaps=do.call(c,listGaps)
Ngaps2=Ngaps[width(Ngaps)>=20]

traceG4 <- list()
traceSH <- list()
for (i in 1:5)
{
Feat <- resize(ori_list_c[[i]],fix='center',width=4001)

pG4 <- G4data[strand(G4data)=='+']
pG4cov <- coverage(pG4)
mG4 <- G4data[strand(G4data)=='-']
mG4cov <- coverage(mG4)

profFeat_p <- pG4cov[Feat]
profFeat_m <- mG4cov[Feat]
prof_G4nearFeat_p <- RleList2matrix(profFeat_p)
prof_G4nearFeat_m <- t(apply(RleList2matrix(profFeat_m),1,rev))
prof_G4nearFeat <- rbind(prof_G4nearFeat_p,prof_G4nearFeat_m)
meantracep <- colMeans(prof_G4nearFeat)
sdtracep <- apply(prof_G4nearFeat,2,sd)
ictopp <- meantracep + (1.96*sdtracep/sqrt(length(Feat)))
icbotp <- meantracep - (1.96*sdtracep/sqrt(length(Feat)))

traceG4[[i]] <- list(prof_G4nearFeat,meantracep,ictopp,icbotp)

altfeat <- shuffleGRgen(1,gen=genome,inputGR=Feat,gap2=Ngaps2,chrlist=1:chnb)
profFeat_p2 <- pG4cov[altfeat]
profFeat_m2 <- mG4cov[altfeat]
prof_G4nearFeat_p2 <- RleList2matrix(profFeat_p2)
prof_G4nearFeat_m2 <- t(apply(RleList2matrix(profFeat_m2),1,rev))
prof_G4nearFeat2 <- rbind(prof_G4nearFeat_p2,prof_G4nearFeat_m2)
meantracep2 <- colMeans(prof_G4nearFeat2)
sdtracep2 <- apply(prof_G4nearFeat2,2,sd)
ictopp2 <- meantracep2 + (1.96*sdtracep2/sqrt(length(Feat)))
icbotp2 <- meantracep2 - (1.96*sdtracep2/sqrt(length(Feat)))

traceSH[[i]] <- list(prof_G4nearFeat2,meantracep2,ictopp2,icbotp2)
}

pdf('./Script_Results/Fig4B.pdf',height=22)
par(mfrow=c(5,1),mai=c(1,0.8,0.4,0.4))
for (i in 1:5)
{
meantrace <- traceG4[[i]][[2]]
ictop <- traceG4[[i]][[3]]
icbot <- traceG4[[i]][[4]]
meantrace2 <- traceSH[[i]][[2]]
ictop2 <- traceSH[[i]][[3]]
icbot2 <- traceSH[[i]][[4]]

xlarg <- 1000L
suff <- 'IS'
ylableg='Fraction of G4FS'
plot(meantrace[(2001L-xlarg):(2001L+xlarg)],x=((-xlarg):xlarg),lwd=2,type='l',xlab=paste('Dist from ',suff,' (bp)',sep=''),ylab='',col='red',ylim=c(0,0.12),xaxs='i',cex.axis=2,cex.lab=2.4,font=2,font.lab=2)
polygon(c(((-xlarg):xlarg),rev((-xlarg):xlarg)),c(ictop[(2001L-xlarg):(2001L+xlarg)],rev(icbot[(2001L-xlarg):(2001L+xlarg)])),border=F,col=rgb(1,0,0,0.2))
lines(meantrace2[(2001L-xlarg):(2001L+xlarg)],x=((-xlarg):xlarg),lwd=1,type='l',col='orange')
polygon(c(((-xlarg):xlarg),rev((-xlarg):xlarg)),c(ictop2[(2001L-xlarg):(2001L+xlarg)],rev(icbot2[(2001L-xlarg):(2001L+xlarg)])),border=F,col=rgb(1,1,0,0.2))
abline(v=0,col='purple',lty=2,lwd=2)
abline(v=-400,col='purple',lty=2,lwd=2)
abline(v=-250,col='purple',lty=2,lwd=2)
legend('topright',c(names(ori_list_c)[i],'shuffled'),text.col=c('red','orange'),bty='n',text.font=2,cex=3)
mtext('Fraction of G4-forming sequence',2,3,cex=1.4,font=2)
}
dev.off()
###

### Fig5B /S6A
# need to annotate mm10 txdb first
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
seqlevels(txdb,pruning.mode='coarse') <- seqlevels(seqinf)
prom <- promoters(txdb,up=2000,down=0)
trans <- transcripts(txdb)
featlist <- list(prom,trans)
allgene <- reduce(do.call(c,featlist),ignore.strand=T)
intergen <- gaps(allgene)
intergen <- intergen[strand(intergen)=='*']
featlist2 <- c(featlist,intergen)
names(featlist2)<- c('Promoter','Transcript','Intergenic')

res3 <- sapply(ori_list_c, function(x) sapply(featlist2, function(y) sum(overlapsAny(x,y,ignore.strand=T))/length(x)))
fFeat <- sapply(featlist2, function(x) sum(width(reduce(x,ignore.strand=T)))/sum(seqlengths(genome)[1:21]))
res4 <- cbind(res3,fFeat)

# There are overlaps between promoters and transcripts so for each class of ori,
# the sum of the ori in the 3 features can be higher than the number of ori in this class

#### Done on a cluster: shuffle ori position and shuffle origroup keeping chromosomal distribution
# DO NOT RUN if not on cluster!
rdlist1000 <- lapply(1:1000,function(x) {shuffleGRgen(x,gen=genome,inputGR=ori_c,gap=Ngaps2,chrlist=1:chnb)})
rdOri_list <- lapply(rdlist1000, function(y) {do.call(c,lapply(1:chnb, function(i) {z<- y[seqnames(y)==seqnames(genome)[i]];z$ori_group <- sample(ori_c[seqnames(ori_c)==seqnames(genome)[i]]$ori_group);return(z)}))})
save(list=c('rdOri_list','rdlist1000'),file='rd_ori.RData')
rdOri_list2 <- sample(rdOri_list,100)
save(rdOri_list2,file='rd_ori2.RData')
# END DO NOT RUN
# A file containing a list of 100 of these shuffled ori is provided
load('rd_ori2.RData')

res.rd <- lapply(rdOri_list, function(z)
{
rd.ori <- split(z,z$ori_group)
res5 <- sapply(rd.ori, function(x) sapply(featlist2, function(y) sum(overlapsAny(x,y,ignore.strand=T))/length(x)))
return(res5)
})
res7 <- array(unlist(res.rd),dim=c(3,5,1000))
res.mean <- apply(res7,c(1,2),mean)
res.sd <- apply(res7,c(1,2),sd)
res6 <- cbind(res.mean,fFeat)

y1 <- apply(res.mean,1,mean)

pdf('./Script_results/Fig5B.pdf')
barplot(rbind(res4[,1],res4[,2],res4[,3],res4[,4],res4[,5],res4[,6]),beside=T,col=c('grey50','red','pink','green','darkgreen','plum'),ylim=c(0,0.7),ylab='Distribution in origin classes',font=2,cex.axis=1.4,cex.lab=1.4,font.lab=2)
legend('topleft',legend=c('insensitive','new','reinforced','reduced','suppressed','genome'),text.col=c('grey50','red','pink','green','darkgreen','plum'),bty='n',cex=1.4,text.font=2)
segments(1,y1[1],7,y1[1],lty=3,lwd=3,col='orange')
segments(8,y1[2],14,y1[2],lty=3,lwd=3,col='orange')
segments(15,y1[3],21,y1[3],lty=3,lwd=3,col='orange')
dev.off()

pdf('./Script_Results/FigS6A.pdf')
bp <- barplot(t(res6),beside=T,col=c('grey50','red','pink','green','darkgreen','plum'),ylim=c(0,0.7),ylab='Distribution in origin classes',font=2,cex.axis=1.4,cex.lab=1.4,font.lab=2)
legend('topleft',legend=c('insensitive','new','reinforced','reduced','suppressed','genome'),text.col=c('grey50','red','pink','green','darkgreen','plum'),bty='n',cex=1.4,text.font=2)
error4barplot(bp,t(res6),t(cbind(res.sd,c(0,0,0))))
dev.off()
###

### Fig7
# data were imported from http://rsat-tagc.univ-mrs.fr/g4/g4_data.html
# done on the genotoul cluster as this operation require large memory
bfname <- 'DMSO_0_virgule_5_pourcent_ACAGTG_L006_R1_PF_sickle_sanger_L30_bowtie2_sorted.bam'
library(GenomicAlignments)
bf <- BamFile(bfname)
ga <- bam2gr(bf)
seqlevels(ga,pruning.mode='coarse') <- seqlevels(G4data)
seqinfo(ga) <- seqinfo(G4data)
rgr <- granges(ga)
rm(ga)
cv <- coverage(rgr)
pG4 <- G4data[strand(G4data)=='+']
mG4 <- G4data[strand(G4data)=='-']

featp <- resize(resize(pG4,fix='start',width=1),fix='center',width=4001)
featm <- resize(resize(mG4,fix='start',width=1),fix='center',width=4001)
profFeat_p <- cv[featp]
profFeat_m <- cv[featm]
prof_G4nearFeat_p <- RleList2matrix(profFeat_p)
meantracep <- colMeans(prof_G4nearFeat_p)
sdtracep <- apply(prof_G4nearFeat_p,2,sd)
ictopp <- meantracep + (1.96*sdtracep/sqrt(dim(prof_G4nearFeat_p)[1]))
icbotp <- meantracep - (1.96*sdtracep/sqrt(dim(prof_G4nearFeat_p)[1]))
rm(prof_G4nearFeat_p)
prof_G4nearFeat_m <- RleList2matrix(profFeat_m)
meantracem <- colMeans(prof_G4nearFeat_m)
sdtracem <- apply(prof_G4nearFeat_m,2,sd)
ictopm <- meantracem + (1.96*sdtracem/sqrt(dim(prof_G4nearFeat_m)[1]))
icbotm <- meantracem - (1.96*sdtracem/sqrt(dim(prof_G4nearFeat_m)[1]))
rm(prof_G4nearFeat_m)
rm(cv)

# Possibility to load precomputed traces
# load("tracereadG4.RData")

pdf('./Script_Results/Fig7.pdf',width=14)
par(mfrow=c(1,2))
xlarg <-2000
ylableg <- 'mean read coverage'
plot(meantracep[(2001L-xlarg):(2001L+xlarg)],x=((-xlarg):xlarg),lwd=2,type='l',xlab="Distance from G4 3' end",col='red',yaxs='i',ylim=c(1,3.2),xaxs='i',cex.axis=1.4,cex.lab=1.8,font=2,font.lab=2,ylab=ylableg)
polygon(c(((-xlarg):xlarg),rev((-xlarg):xlarg)),c(ictopp[(2001L-xlarg):(2001L+xlarg)],rev(icbotp[(2001L-xlarg):(2001L+xlarg)])),border=F,col=rgb(1,0,0,0.2))
abline (v=0,lty=5,lwd=4,col='grey',lend=1)
text(-350,3,'G4',font=2,cex=2.5)
plot(meantracem[(2001L-xlarg):(2001L+xlarg)],x=((-xlarg):xlarg),lwd=2,type='l',xlab="Distance from G4 3' end",col='blue',yaxs='i',ylim=c(1,3.2),xaxs='i',cex.axis=1.4,cex.lab=1.8,font=2,font.lab=2,ylab=ylableg)
polygon(c(((-xlarg):xlarg),rev((-xlarg):xlarg)),c(ictopm[(2001L-xlarg):(2001L+xlarg)],rev(icbotm[(2001L-xlarg):(2001L+xlarg)])),border=F,col=rgb(0,0,1,0.2))
abline (v=0,lty=5,lwd=4,col='grey',lend=2)
text(350,3,'G4',font=2,cex=2.5)

dev.off()
###

### FigS4D
pdf('./Script_Results/FigS4D.pdf')
boxplot(list(ori_list_c[[1]]$ori_FC,ori_list_c[[2]]$ori_FC,ori_list_c[[3]]$ori_FC,ori_list_c[[4]]$ori_FC,ori_list_c[[5]]$ori_FC),outline=F,ylab='',names=rep("",5),cex.axis=1.6)
mtext('Log2(Fold Change)',2,2.5,font=2,cex=1.8)
mtext(c('Insensitive','New','Reinforced','Reduced','Suppressed'),1,1,at=1:5,font=2,cex=1.2)
axis(2,font.axis=2,cex.axis=1.6)
dev.off()

###

### FigS5B
G4_1kb_ori_c <- endoapply(ori_list_c, function(x) G4data[overlapsAny(G4data,x+500)])

mxscore <- lapply(G4_1kb_ori_c, function(x) abs(x$score))
pdf('./Script_Results/FigS5B.pdf')
par(mfrow=c(2,2))
b <- c(25,25,12,25)
dum <-lapply(c(1,2,4,3),function(i) hist(mxscore[[i]], main=toupper(names(mxscore)[i]),xlim=c(1,4),xlab='G4H_score',breaks=b[i],font=2,font.lab=2,font.axis=2,cex.main=1.4,cex.lab=1.4,cex.axis=1.2))
dev.off()
###

### FigS5C
pdf('./Script_Results/FigS5C.pdf')
boxplot(lapply(G4_1kb_ori_c, width),outline=F,ylim=c(0,100),ylab='',names=rep("",5),cex.axis=1.6)
mtext('G4 width (bp)',2,2.5,font=2,cex=1.8)
mtext(c('Insensitive','New','Reinforced','Reduced','Suppressed'),1,1,at=1:5,font=2,cex=1.2)
axis(2,font.axis=2,cex.axis=1.6)
dev.off()
###

### FigS5D
ori_reads_G4h2 <- flankingG4H(ori_reads,G4data,1000)
test <- ori_reads_G4h2[ori_reads_G4h2$ori_group!='new']
test_red <- test
test_red$G4class_preIS <- test$preIS
test_red$G4class_preIS[test_red$G4class_preIS>5] <- 6
test_red$G4class_postIS <- test$postIS
test_red$G4class_postIS[test_red$G4class_postIS>5] <- 6

pdf('./Script_Results/FigS5D.pdf',width=12)

bp1 <- boxplot(test_red$ori_DMSreads~test_red$G4class_preIS,outline=F,main=paste0('All origins tested - ',length(test)),xlab='',ylab='',ylim=c(0,800),names=rep('',7),cex.axis=1.4,cex.main=1.8)
mtext('Origin strength (Reads number)',2,2.5,font=2,cex=1.8)
mtext(c('0','1','2','3','4','5','6+'),1,1,at=1:7,font=2,cex=1.4)
axis(2,font.axis=2,cex.axis=1.4)
mtext (paste0('n = ',bp1$n),1,2,at=1:7,font=4,cex=1)
mtext('Number of G4s upstream to the IS',1,3.5,font=2,cex=1.8)

# alternative version
bp1 <- boxplot(test_red$ori_DMSreads~test_red$G4class_preIS,outline=F,main=paste0('All origins tested n = ',length(test)),xlab='',ylab='',ylim=c(0,800),names=rep('',7),cex.axis=1.4,cex.main=1.8)
mtext('Origin strength (Reads number)',2,2.5,font=2,cex=1.8)
mtext(c('0','1','2','3','4','5','6+'),1,1,at=1:7,font=2,cex=1.4)
axis(2,font.axis=2,cex.axis=1.4)
mtext (paste0('n = ',bp1$n),1,2,at=1:7,font=4,cex=1,col='blue')
mtext('Number of G4s upstream to the IS',1,3.5,font=2,cex=1.8)
par(new=T)
barplot((bp1$n)/25, col=rgb(0,0,1,0.3),axes=F,add=T,width=0.2,space=4,border=rgb(0,0,1,0.1))

dev.off()
###


### TableS1
# G4H score have been computed using the G4HunterTable shiny app reported in Lacroix et al Bioinformatics 2018.
# https://github.com/LacroixLaurent/G4HunterTable
###

### TableS3
# G4H score have been computed using the G4HunterTable shiny app reported in Lacroix et al Bioinformatics 2018.
# https://github.com/LacroixLaurent/G4HunterTable
###

### TableS4
# beware if your setup is not multicore friendly
# QP37 <- do.call(c,mclapply(1:chnb, quadparser))
# Possibility to load precomputed quadparser result with default parameters
load("mm10_QP37.RData")
QP37 <- GRmm10_QP37
ori <- ori_reads
ori$QP <- overlapsAny(resize(ori,fix='center',width=1001),QP37)
ori$G4H2 <- overlapsAny(resize(ori,fix='center',width=1001),G4data)

ori3 <- flankingG4H(ori_reads,GRmm10_QP37,501)
res5 <- cbind(table(ori3$ori_group),table(ori3[ori3$preIS!=0]$ori_group),table(ori3$ori_group)-table(ori3[ori3$preIS!=0]$ori_group))

ori2 <- flankingG4H(ori_reads,G4data,501)
res3 <- cbind(table(ori2$ori_group),table(ori2[ori2$preIS!=0]$ori_group),table(ori2$ori_group)-table(ori2[ori2$preIS!=0]$ori_group))

res.table <- cbind(res5,res3[,2:3])
res.table2 <- rbind(res.table, colSums(res.table),colSums(res.table[c(1,3:5),]),colSums(res.table[c(1:4),]),colSums(res.table[c(2:5),]))
res.table3 <- cbind(res.table2,round(res.table2[,2]/res.table2[,1]*100,1),round(res.table2[,4]/res.table2[,1]*100,1))
colnames(res.table3) <- c('total','QP+','QP-','G4H+','G4H-','%QP+','%G4H+')
rownames(res.table3)[6:9] <- c('total','all control','all PhenDC3','oris changed')
write.table(res.table3,file='./Script_Results/TableS4.txt', sep='\t',col.names=NA)
###

### Fig4D
pdf('./Script_Results/Fig4D.pdf',width=12)
par(mai=c(0.4,1,0.4,0.4))
barplot(rbind(res.table3[c(6,1:5),7],100-res.table3[c(6,1:5),7]),beside=F,col=c('grey42','grey66'),ylab='', names=rep('',6),border='white',las=1,ylim=c(0,115),font.axis=2,cex.axis=1.4,space=1)
legend('top', fill=c('grey42','grey66'),legend=c('G4+','G4-'),bty='n',cex=2,horiz=T,text.font=2,border='white')
mtext('Distribution in origin classes (%)',2,3,font=2,cex=2)
mtext(c('Total','Insensitive','New','Reinforced','Reduced','Suppressed'),1,1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),font=2,cex=1.8)

dev.off()
###
