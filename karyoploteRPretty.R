library("karyoploteR")
library("regioneR")
library("GenomicAlignments")
library("zoo")
library("argparse")
#--------------------------------------------------------
parser <- ArgumentParser(description='Process some integers')
#--------------------------------------------------------
parser$add_argument('--chr', 
                    metavar='N', 
                    type="character", 
                    nargs='+',
                    default="chr9",
                    help='a chromosome id')
#--------------------------------------------------------
parser$add_argument('--file', 
                    dest='infile', 
                    help='input files (default: find the max)')
#--------------------------------------------------------
parser$add_argument('--tbam', 
                    dest='tbam',
                    default='max',
                    help='tumor bam to store coverage')
#--------------------------------------------------------
parser$add_argument('--nbam', 
                    dest='infile', 
                    default='max',
                    help='normal bam to store coverage')
#--------------------------------------------------------
parser$print_help()
# default args for ArgumentParser()$parse_args are commandArgs(TRUE)
args <- parser$parse_args()
infile<-get(args$infile)
gmax  <-get(args$max)
chrs  <-get(args$chr)
# default args for ArgumentParser()$parse_args are commandArgs(TRUE)
tbam <- get(args$tbam)
nbam <- get(args$nbam)
#--------------------------------------------------------
tcvg <- coverage(readGAlignments(tbam))
ncvg <- coverage(readGAlignments(nbam))
#--------------------------------------------------------
d<-read.table("H:/projects/ctDNA_PROJECTS/test.txt",
              header=T)
gmax<-2
#--------------------------------------------------------
dp.colors  <- c("#FFBD07", "#00A6ED",  "#FF3366", "#8EA604", "#C200FB")
data.points<- GRanges(IRanges(start=d$start,end=d$end),
                      seqnames=d$chromosome,
                      gene=d$gene,
                      y=d$log2,
                      pvalue=d$pvalue)
#--------------------------------------------------------
data.ranges=data.points
data.points.colors=dp.colors
num.data.points      <- 3000
num.big.regions.up   <- 30
num.big.regions.down <- 30
num.mid.regions      <- 6000
num.marks            <- 90
#--------------------------------------------------------
gains <-d[which(d$log2>0),]
losses<-d[which(d$log2<0),]
#--------------------------------------------------------
data.ranges<-GRanges(chr=d$chr,
                     IRanges(start=d$start,end=d$end),
                     seqnames=d$chr,genes=d$gene,
                     logR=d$log2,pvalue=d$pvalue)
#--------------------------------------------------------
gains.genes     <-unique(gains$gene[which(gains$fdr<0.1)])
gains.genes.data<-data.ranges[match(gains.genes,d$gene)]
loss.genes      <-unique(losses$gene[which(losses$fdr<0.1)])
loss.genes.data <-data.ranges[match(loss.genes,d$gene)]
#--------------------------------------------------------
mid.regs <- toGRanges(d[which(d$pvalue >= 0.05),])
#--------------------------------------------------------
kp<-plotKaryotype(genome="hg38",
                  chromosomes=c(chrs),
                  plot.type=2)
#--------------------------------------------------------
### Data Panel 1 ###
#--------------------------------------------------------
#   Big regions
#--------------------------------------------------------
big.regs.up<-gains.genes.data
big.regs.down<-loss.genes.data
#--------------------------------------------------------
kpRect(kp, data = big.regs.up, 
       y0=0, y1=1, col="#FFDDDD", 
       border=NA, r0=0, r1=0.8)
#--------------------------------------------------------
kpRect(kp, data = big.regs.down, 
       y0=0, y1=1, col="#DDFFDD", 
       border=NA, r0=0, r1=0.8)
#--------------------------------------------------------
#Data points
yrange<-max(data.points$y)-min(data.points$y)
#--------------------------------------------------------
kpAxis(kp, 
       ymin = -25, 
       ymax = 5, 
       r0=0.2, r1=1, numticks = 7, 
       col="#666666", cex=0.5)
#--------------------------------------------------------
yh<-1/yrange*abs(min(data.points$y))
#--------------------------------------------------------

#Mean and sd of the data points.  
for(chr in seqlevels(kp$genome)) 
{
  chr.dp <- sort(keepSeqlevels(x = data.points, 
                               value = chr, 
                               pruning.mode = "coarse"))
  
  rmean <- zoo::rollmean(chr.dp$y, 
                         k = 5, align = "center")  
  
  rsd <- zoo::rollapply(data = chr.dp$y, 
                        FUN=sd, width=8)
  
  kpLines(kp, 
          chr = chr, 
          x=start(chr.dp)[3:(length(chr.dp)-3)], 
          y=rmean,
          ymin=min(rmean),
          ymax=max(rmean),
          col=data.points.colors[3])
  
  kpPlotRibbon(kp, chr=chr, data=chr.dp[3:(length(chr.dp)-3)], 
               y0=rmean-rsd, 
               y1=rmean+rsd, 
               ymin=min(rmean-rsd),
               ymax=max(rmean+rsd),
               col="#FF336633",
               border=NA)
}

#--------------------------------------------------------

signif.points=data.points[which(data.points$pvalue<0.05),]

kpPoints(kp, 
         data=signif.points,
         y=signif.points$y,
         pch=16, cex=0.5, 
         col=rgb(0.5,0,0,0.4), 
         r0=0.5, r1=0.8)

kpAddLabels(kp, 
            labels = "logR", 
            srt=90, pos=1, 
            label.margin = 0.04, 
            ymax=ymax, 
            ymin=ymin)

signif.points<-signif.points[seqnames(signif.points)==chr,]
signif.points<-signif.points[order(signif.points$pvalue),]
top.genes    <-reduce(split(signif.points, 
                            elementMetadata(signif.points)$gene),
                      drop.empty.ranges=TRUE)
top.genes    <-top.genes[which(lapply(top.genes,length)>0)]

top.pvals    <-lapply(split(signif.points,
                            signif.points$gene),function(x){return(min(x$pvalue,na.rm=T))})

top.pvals    <-order(unlist(top.pvals[match(names(top.genes),names(top.pvals))]))[1:params$N]

gene.names<-names(top.genes)[top.pvals]
top.genes<-top.genes[top.pvals]

kpPlotMarkers(kp,
              chr=chr,
              x=unlist(lapply(top.genes,function(x){return(min(ranges(x)@start))})),
              labels = gene.names, 
              r0=0.9,
              adjust.label.position=T,
              label.margin=0.01,
              text.orientation = "horizontal")

### Data Panel 2 ###
# bam coverage plot data ###
#medium regions and their coverage
kpPlotRegions(kp,  
              data=cvg, 
              r0 = 0.2, r1=1, 
              border=NA, 
              data.panel=2)

kpPlotCoverage(kp, 
               data=cvg, 
               r0=0.2, r1=0, 
               col=data.points.colors[2], 
               data.panel = 2)

kpPlotCoverage(kp, 
               data=cvg, 
               r0=0.2, r1=0.12, 
               col=data.points.colors[1], 
               data.panel = 2)

kpText(kp, chr=seqlevels(kp$genome), 
       y=0.4, x=0, data.panel = 2, r0=0.2, r1=0, 
       col="#444444", label="30x", cex=0.8, pos=2)

kpAbline(kp, h=0.4, data.panel = 2, 
         r0=0.2, r1=0, col=data.points.colors[3])



#---------------------------------------------------
# Data Panel 3 #
# Highlight top N genes ###
#---------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
txdb<-genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
map<-as.list(org.Hs.egSYMBOL)
m<-match(gene.names,map)
txdb<-txdb[match(names(map[m]),names(txdb))]
#---------------------------------------------------
# top N genes panel
# mark by copy number???
#---------------------------------------------------
kpPlotGenes(kp, 
            data=txdb, 
            plot.transcripts.structure=FALSE, 
            add.transcript.names=FALSE, 
            gene.names=gene.names, 
            r1=0.8, 
            col="blue", 
            marks.col="white", 
            gene.name.col="black")

#--- intersect genomic intervals ----


