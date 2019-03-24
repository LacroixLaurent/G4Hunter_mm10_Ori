################################################################################
### Function related to G4Hunter_mm10_Ori
### by L. Lacroix, laurent.lacroix@inserm.fr
### 20190323
################################################################################


################################################################################
######################### G4Hunter functions ###################################
################################################################################

###### G4translate change the DNA code into G4Hunter code.
###### Only G or C are taken into account. non G/C bases are translated in 0
###### It is OK if N or U are in the sequence, but might be a problem if other letters or numbers are present
G4translate <- function(y,v1=1,v2=2,v3=3,v4=4)		
# x a DNAString or a DNAStringSet (or just a string of char)
{
require(XVector)
x= toupper(Rle(strsplit(as.character(y),NULL)[[1]]))
xres=x
runValue(xres)[runValue(x)=='C' & runLength(x)>3] <- -v4
runValue(xres)[runValue(x)=='C' & runLength(x)==3] <- -v3
runValue(xres)[runValue(x)=='C' & runLength(x)==2] <- -v2
runValue(xres)[runValue(x)=='C' & runLength(x)==1] <- -v1
runValue(xres)[runValue(x)=='G' & runLength(x)>3] <- v4
runValue(xres)[runValue(x)=='G' & runLength(x)==3] <- v3
runValue(xres)[runValue(x)=='G' & runLength(x)==2] <- v2
runValue(xres)[runValue(x)=='G' & runLength(x)==1] <- v1
runValue(xres)[runValue(x)!='C' & runValue(x)!='G'] <- 0
# N or U are not a problem
Rle(as.numeric(xres))
}

##### return the G4Hscore. y is a string of char (DNA sequence)
G4Hscore <- function(y)
	{
		y3=G4translate(y)
		mean(y3)
	}

##### functions to refine G4hunt results
G4startrun=function(y,chrom=chr,letter='C')	#y is a START
	{
		if (letter(chrom,y)==letter)
			{
			while (letter(chrom,y)==letter & y!=1) 
				{
					if (letter(chrom,y-1)==letter) 
					{
						y <- y-1
					}else{
						break
					}
				}
			}else{
			y=y+1
			while (letter(chrom,y)!=letter) {y=y+1}
			}	
		return(y)
	}

G4endrun=function(y,chrom=chr,letter='C')	#y is a END
	{
		if (letter(chrom,y)==letter)
			{
			while (letter(chrom,y)==letter & y!=length(chrom))
				{
					if (letter(chrom,y+1)==letter)
					{
						y <- y+1
					}else{
						break
					}
				}
			}else{
			y=y-1
			while (letter(chrom,y)!=letter) {y=y-1}
			}
		return(y)
	}

#### G4Hunt function with refining and max score
mG4huntref <- function(i,k=25,hl=1.5,gen=genome,masked=5,with.seq=T,Gseq.only=T)
{
require(GenomicRanges,Biostrings)
# k=RUNMEAN WINDOW SIZE, hl=threshold, i=chromosom number in the genome

chr=gen[[i]]
if (masked==2) {chr=injectHardMask(chr)}
if (masked==3) {active(masks(chr))['RM']=T;chr=injectHardMask(chr)}
if (masked==0) {active(masks(chr))=F;chr=injectHardMask(chr)}
if (masked==4) {active(masks(chr))=T;chr=injectHardMask(chr)}

tchr <- G4translate(chr)
chr_G4hk <- runmean(tchr,k)
if (class(gen)=="DNAStringSet")
	{seqname=names(gen)[i]
	}else{
	seqname=seqnames(gen)[i]
	}	

j <- hl
chrCh <- Views(chr_G4hk, chr_G4hk<=(-j))
chrGh <- Views(chr_G4hk, chr_G4hk>=j)

IRC <- reduce(IRanges(start=start(chrCh),end=(end(chrCh)+k-1)))
if (length(IRC)==0)
	{
	nxC <- GRanges()                                
	}else{
	nnIRC=IRC
	start(nnIRC)=sapply(start(IRC),G4startrun,letter='C',chrom=chr)
	end(nnIRC)=sapply(end(IRC),G4endrun,letter='C',chrom=chr)
	seqC=as.character(Views(chr,nnIRC))
	if (Gseq.only) 
		{
			nnseqC=as.character(reverseComplement(Views(chr,nnIRC)))
		}else{
			nnseqC=as.character(Views(chr,nnIRC))
		}
	nG4scoreC=sapply(seqC,function(x) signif(G4Hscore(x),3))
	mscoreC <- signif(min(Views(chr_G4hk,IRC)),3)
	straC <- Rle(rep('-',length(IRC)))
	hlC <- Rle(rep(j,length(IRC)))
	kC <- Rle(rep(k,length(IRC)))
	maskC <- Rle(rep(masked,length(IRC)))      
	nxC <- GRanges(seqnames=Rle(seqname),  
    ranges=nnIRC,
    strand=straC,
    score=nG4scoreC,
    max_score=mscoreC,
    hl=hlC,
    k=kC,
    mask=maskC,
    sequence=nnseqC,
    seqinfo=seqinfo(gen))
	}

IRG <- reduce(IRanges(start=start(chrGh),end=(end(chrGh)+k-1)))
if (length(IRG)==0)
	{
	nxG <- GRanges()                                
	}else{
	nnIRG=IRG
	start(nnIRG)=sapply(start(IRG),G4startrun,letter='G',chrom=chr)
	end(nnIRG)=sapply(end(IRG),G4endrun,letter='G',chrom=chr)
	nnseqG=as.character(Views(chr,nnIRG))
	nG4scoreG=sapply(nnseqG,function(x) signif(G4Hscore(x),3))
	mscoreG <- signif(max(Views(chr_G4hk,IRG)),3)
	straG <- Rle(rep('+',length(IRG)))
	hlG <- Rle(rep(j,length(IRG)))
	kG <- Rle(rep(k,length(IRG)))
	maskG <- Rle(rep(masked,length(IRG)))      
	nxG <- GRanges(seqnames=Rle(seqname),  
    ranges=nnIRG,
    strand=straG,
    score=nG4scoreG,
    max_score=mscoreG,
    hl=hlG,
    k=kG,
    mask=maskG,
    sequence=nnseqG,
    seqinfo=seqinfo(gen))
	}

nx <- sort(c(nxC,nxG),ignore.strand=T)
names(nx) <- NULL
if (with.seq==F) {nx$sequence=NULL}
return(nx)
}

### my version of quadparser
quadparser <- function(i,gen=genome,gmin=3,gmax='',lmin=1,lmax=7)
# this function extract form a chromosome 'i' in a genome from a BSGenome all the sequence matching the Gpattern
# "G{gmin,gmax}.{lmin,lmax}G{gmin,gmax}.{lmin,lmax}Ggmin,gmax}.{lmin,lmax}G{gmin,gmax}"
# or its complement pattern
# by default, the pattern is "G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}"
# the * option allow to compensate for overlapping excluded from gregexpr
	{
		require(GenomicRanges)

		Gpat <- paste('G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}.{',lmin,',',lmax,'}G{',gmin,',',gmax,'}','(.{',lmin,',',lmax,'}G{',gmin,',',gmax,'})*',sep='')
		Cpat <- paste('C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}.{',lmin,',',lmax,'}C{',gmin,',',gmax,'}','(.{',lmin,',',lmax,'}C{',gmin,',',gmax,'})*',sep='')


	if (class(gen)=="DNAStringSet")
	{seqname=names(gen)[i]
	}else{
	seqname=seqnames(gen)[i]
	}	
		Gres=gregexpr(Gpat, gen[[i]], perl=T)
		Cres=gregexpr(Cpat, gen[[i]], perl=T)
		
		if (Gres[[1]][1]>0) 
			{Gwidth=attr(Gres[[1]],"match.length")
			Gstart=unlist(Gres)
			G4Ranges=GRanges(seqnames=seqname,IRanges(start=Gstart,width=Gwidth),strand='+',seqinfo=seqinfo(gen))
			} else {
			G4Ranges=GRanges()}
		if (Cres[[1]][1]>0)		
			{Cwidth=attr(Cres[[1]],"match.length")
			Cstart=unlist(Cres)
			C4Ranges=GRanges(seqnames=seqname,IRanges(start=Cstart,width=Cwidth),strand='-',seqinfo=seqinfo(gen))
			}else{
			C4Ranges=GRanges()}	
		res <- sort(c(G4Ranges,C4Ranges),ignore.strand=T)
		res$seq <- as.character(Views(gen[[i]],ranges(res)))
		return(res)
	}


################################################################################

################################################################################
### helper from PGP Martin 
RleList2matrix <- function(rlelist)
# rlelist should be a RleList, not just a list of Rle
{
matrix(as.numeric(unlist(rlelist,use.names=F)),
        nrow=length(rlelist),
        byrow=T,
        dimnames=list(names(rlelist),NULL))
}
################################################################################

################################################################################
### other functions by LLacroix

# Importing data
# for SNS peak data
rawread2gr <- function(x,seqinfo=seqinf)
{
require('GenomicRanges')

x=GRanges(seqnames=Rle(x[,1]),
        ranges=IRanges(start=x[,2],end=x[,3]),
        strand=Rle(rep("*",nrow(x))),
        ori_group=x[,15],seqinfo=seqinfo,ori_FC=x[,14],ori_DMSreads=x[,11],ori_PhDCreads=x[,12])
return(x)
}

# for bam file
bam2gr <- function(bf,paired=F)
{
require(GenomicAlignments)
		if (paired) 
		{ 
			scanpar2 <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE,isProperPair=TRUE))
			# keep only properly paired reads
			ga <- readGAlignmentPairs(bf,param=scanpar2)
			gra <- granges(ga)
			gra <- gra[!duplicated(gra)]
			# remove duplicated fragments
			ga <- gra
		}
		else
		{
			scanpar <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE))
			ga <- readGAlignments(bf,param=scanpar)
			dup <- duplicated(granges(ga))
			ga <- ga[!dup]
			# remove duplicated reads
		}
	return(ga)
}

################################################################################
# looking for G4 near ori
flankingG4H <- function (oriGR, featdata=G4data2,distance=1000)
{
	toto <- resize(oriGR, fix='center', width=1)
	featdata_p <- featdata[strand(featdata)=='+']
	featdata_m <- featdata[strand(featdata)=='-']
	
	oriGR$pre_p <- countOverlaps(resize(toto,fix='end',width=distance),featdata_p)
	oriGR$pre_m <- countOverlaps(resize(toto,fix='end',width=distance),featdata_m)
	oriGR$fol_p <- countOverlaps(resize(toto,fix='start',width=distance),featdata_p)
	oriGR$fol_m <- countOverlaps(resize(toto,fix='start',width=distance),featdata_m)
	oriGR$preIS <- oriGR$pre_p+oriGR$fol_m
	oriGR$postIS <- oriGR$pre_m+oriGR$fol_p
	oriGR$center <- countOverlaps(toto,featdata)
	return(oriGR)
}

################################################################################


################################################################################
### GenomicRanges Shuffling function
################################################################################

# function to resample on a given chromosom keeping strand information
shuffleGR3=function(gen=genome,chrnb=24,inputGR=G4data,gap=Ngaps2)
	{	require(GenomicRanges)
		hit <- inputGR[seqnames(inputGR)==seqnames(gen)[chrnb]]
		gapchr=gap[seqnames(gap)==seqnames(gen)[chrnb]]
		rgap <- ranges(gapchr)
		ravail <- gaps(rgap)
#		st_avail <- unlist(as.vector(ravail))
# broken in BioC3.7, should come back in BioC3.8
# Temporary fix
		st_avail <- IRanges:::unlist_as_integer(ravail)
#		
		st_rdgr <- sample(st_avail,length(hit))
		if (length(hit)==1)		
				{
				wi_rdgr <- width(hit)
				}else{
				wi_rdgr <- sample(width(hit))
				#necessary if only one range sample(width()) choose a number
				#betwen in 1:width() rather than one width
				}
		ra_rdgr <- sort(IRanges(start=st_rdgr,width=wi_rdgr))
		keep <- IRanges()
		ra_rdgr2 <- IRanges()
		while ((sum(overlapsAny(ra_rdgr,rgap))!=0) | (sum(overlapsAny(ra_rdgr2,keep))!=0))
			{
			keep <- ra_rdgr[overlapsAny(ra_rdgr,rgap)==0]
			hit2 <- ra_rdgr[overlapsAny(ra_rdgr,rgap)!=0]
			st_rdgr2 <- sample(st_avail,length(hit2))
			if (length(hit2)==1)
				{
				wi_rdgr2 <- width(hit2)
				}else{
				wi_rdgr2 <- sample(width(hit2))
				}
			ra_rdgr2 <- IRanges(start=st_rdgr2,width=wi_rdgr2)	
			ra_rdgr <- c(keep,ra_rdgr2)
			}
		rdgr <- sort(GRanges(seqnames=Rle(rep(seqnames(gen)[chrnb],length(hit))),ranges=ra_rdgr,strand=Rle(rep('*',length(hit))),seqinfo=seqinfo(inputGR)))
		return(rdgr)
	}

# function to resample on a genome

shuffleGRgen <- function(dummy=1,gen2=genome,inputGR2=G4data,gap2=Ngaps2,chrlist=1:chnb)
	{
	rdlist=GRangesList()
	for (i in chrlist) {rdlist[[i]] <- shuffleGR3(gen=gen2,chrnb=i,inputGR=inputGR2,gap=gap2)}
	y<- do.call(c,rdlist)
	return(y)
	}

# Gap annotation
findNgaps <- function(x)
# x is a DNAString
	{ y=Rle(strsplit(as.character(x),NULL)[[1]])
	  y2=ranges(Views(y,y=='N'))
	  return(y2)	# y2 is a list of IRanges
	}

################################################################################

################################################################################