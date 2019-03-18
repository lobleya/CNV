#-----------------------------------------------------------------
library("karyoploteR")
library("regioneR")
library("zoo")
#library("GenomicAlignments")
library("argparse")
#-----------------------------------------------------------------
parser <- ArgumentParser(description='make pretty karyotype plots')
parser$add_argument('--chr', 
                    metavar='chr', 
                    type="character", 
                    nargs='+',
                    default="chr9",
                    help='select chromosome for plotting')
#--------------------------------------------------------
#parser$add_argument('--tbam', dest='tbam', type="character",
#                    help='bamfile for coverage plots',
#                    default="H://projects/ctDNA_pancreas/tbam")
#--------------------------------------------------------
#parser$add_argument('--nbam', dest='nbam', type="character",
#                    help='bamfile for coverage plots',
#                    default="H://projects/ctDNA_pancreas/nbam")
#--------------------------------------------------------
parser$add_argument('--logR', dest='logRfile', type="character",
                    help="logR data file for processing",
                    default="H://projects/ctDNA_pancreas/test.txt")
#--------------------------------------------------------
parser$add_argument('--topN', dest="N",type="integer",
                    help="top N genes for display", default=10)
#--------------------------------------------------------
parser$print_help()
N<-parser$N
#--------------------------------------------------------
# default args for ArgumentParser()$parse_args are
#--------------------------------------------------------
parser$logR<-"H://projects/ctDNA_PROJECTS/test.txt"
#--------------------------------------------------------
d<-read.table(parser$logR,
              header=T)
d<-d[grep("Antitarget",d$gene,invert=T),]
d$pvalue<-1-pnorm(abs(d$log2-mean(d$log2))/sd(d$log2))
#--------------------------------------------------------
gmax<-2
chrs<-parser$chr
#--------------------------------------------------------
#tcvg <- coverage(readGAlignments(parser$tbam))
#ncvg <- coverage(readGAlignments(parser$nbam))
#--------------------------------------------------------
dp.colors  <- c("#FFBD07", "#00A6ED",  "#FF3366", "#8EA604", "#C200FB")
#-----------------------------------------------------------------
data.points<- GRanges(IRanges(start=d$start,end=d$end),
                      seqnames=d$chromosome,
                      gene=d$gene,
                      y=d$log2,
                      pvalue=d$pvalue)

#--------------------------------------------------------
signif.points<-toGRanges(d[which(d$pvalue<0.05 & d$chromosome==chr),])
smin<-min(signif.points$log2)
smax<-max(signif.points$log2)
#-----------------------------------------------------------------
data.ranges=data.points
data.points.colors=dp.colors
num.data.points      <- 3000
num.big.regions.up   <- 30
num.big.regions.down <- 30
num.mid.regions      <- 6000
num.marks            <- 90
#-----------------------------------------------------------------
pdf(paste(chr,".pdf",sep=""),width=10,height=6)
#-----------------------------------------------------------------
gains <-d[which(d$log2>0),]
losses<-d[which(d$log2<0),]
#-----------------------------------------------------------------
data.ranges<-GRanges(chr=d$chr,
                     IRanges(start=d$start,end=d$end),
                     seqnames=d$chr,genes=d$gene,
                     logR=d$log2,pvalue=d$pvalue)
#--------------------------------------------------------
gains.genes<-unique(gains$gene[which(gains$fdr<0.1)])
gains.genes.data<-data.ranges[match(gains.genes,d$gene)]
loss.genes<-unique(losses$gene[which(losses$fdr<0.1)])
loss.genes.data<-data.ranges[match(loss.genes,d$gene)]
#--------------------------------------------------------
mid.regs <- GRanges(d[which(d$pvalue >= 0.05 & d$chr==chr),])
#--------------------------------------------------------
kp<-plotKaryotype(genome="hg38",
                  chromosomes=chr,
                  plot.type=2)
#-----------------------------------------------------------------
### Data Panel 1 ###
#-----------------------------------------------------------------
#Big regions
#-----------------------------------------------------------------
big.regs.up<-GRanges(d[which(d$pvalue <0.05 & d$chr==chr & d$log2>0 ),])
big.regs.up$x0<-big.regs.up@ranges@start
big.regs.up$x1<-big.regs.up$x0+big.regs.up@ranges@width
#-----------------------------------------------------------------
big.regs.down<-GRanges(d[which(d$pvalue <0.05 & d$chr==chr & d$log2 < 0),])
big.regs.down$x0<-big.regs.down@ranges@start
big.regs.down$x1<-big.regs.down$x0+big.regs.down@ranges@width
#--------------------------------------------------------
kpRect(kp, data = mid.regs, 
        col="lightgray", 
        y0=0.5,
        y1=0.7,
       border="lightgray", 
       r0=0, 
       r1=1,
       data.panel = 1)
#--------------------------------------------------------
#Data points
yrange<-max(data.points$y)-min(data.points$y)
#--------------------------------------------------------
ymin<-  25
ymax<- -25
#--------------------------------------------------------
kpAxis(kp, 
       ymin = smin, 
       ymax = smax, 
       r0=0.2, 
       r1=1, 
       numticks = 11, 
       col="#666666", cex=0.5)
#--------------------------------------------------------
yh<-1/yrange*abs(min(data.points$y))
#--------------------------------------------------------
#Mean and sd of the data points.  
chr.dp <- sort(keepSeqlevels(x = data.points, 
                             value = chr, 
                             pruning.mode = "coarse"))

rmean <- zoo::rollmean(chr.dp$y, 
                       k = 5, align = "center")  
  
rsd <- zoo::rollapply(data = chr.dp$y, 
                      FUN=sd, width=8)
#--------------------------------------------------------

kpLines(kp, 
        chr =chr, 
        x   =start(chr.dp)[3:(length(chr.dp)-3)], 
        y   =rmean,
        ymin=smin,
        ymax=smax,
        col =dp.colors[3], 
        r0  =0.3, 
        r1  =1)
#--------------------------------------------------------
kpPlotRibbon(kp, 
             chr=chr, 
             data=chr.dp[3:(length(chr.dp)-3)],
             y0=rmean-rsd,
             y1=rmean+rsd,
             ymin=smin,
             ymax=smax,
             r0=0.2, 
             col="#FF336633", 
             border=NA)
#--------------------------------------------------------
kpPoints(kp, 
         data=signif.points,
         y=signif.points$log2*0.9,
         pch=16, 
         cex=0.5, 
         col="brown", 
         r0=0.3, 
         r1=1, 
         ymin=smin,
         ymax=smax)
#--------------------------------------------------------
ymax=max(rmean)
ymin=min(rmean)
#--------------------------------------------------------
kpAddLabels(kp, col="red",
            labels = "logR", 
            pos=1, 
            label.margin = 0.15,
            offset=0, 
            r0=ymax, 
            r1=ymin,
            data.panel=1)
#--------------------------------------------------------
signif.points<-signif.points[order(signif.points$pvalue),]
signif.genes<-signif.points[grep("-",signif.points$gene,invert=T),]
up.signif.genes<-signif.genes[which(signif.genes$log2>0),]
down.signif.genes<-signif.genes[which(signif.genes$log2<0),]
top.genes    <-reduce(split(signif.genes, 
                            elementMetadata(signif.genes)$gene),
                            drop.empty.ranges=TRUE)
top.genes    <-top.genes[which(lapply(top.genes,length)>0)]
                      
top.pvals    <-lapply(split(signif.genes,
                               signif.genes$gene),function(x){return(min(x$pvalue,na.rm=T))})

top.pvals    <-order(unlist(top.pvals[match(names(top.genes),
                                            names(top.pvals))]))[1:N]
#---------------------------------------------------------
gene.names<-names(top.genes)[top.pvals]
top.genes <-top.genes[top.pvals]
up.genes  <-top.genes[which(names(top.genes) %in% up.signif.genes$gene), ]
down.genes<-top.genes[which(names(top.genes) %in% down.signif.genes$gene), ]
#--------------------------------------------------------
if(length(up.genes)>1)
{
kpPlotMarkers(kp,
              chr=chr,
              x=unlist(lapply(up.genes,function(x){return(min(ranges(x)@start))})),
              labels = names(up.genes), 
              r0=0.8,
              r1=1.2,
              offset=0.5,
              line.color="gray",
		  label.color="gray",
              adjust.label.position=T,
              label.margin=0.8,
              text.orientation = "horizontal",
              cex=0.6)
}

if(length(down.genes)>1)
{
kpPlotMarkers(kp,
              chr=chr,
              x=unlist(lapply(down.genes,function(x){return(min(ranges(x)@start))})),
              labels = names(down.genes), 
              r0=0.4,
              r1=0,
              offset=-0.5,
		          label.color="gray",
              line.color="gray",
              adjust.label.position=T,
              label.margin=0.8,
              text.orientation = "horizontal",
              cex=0.6)
}
#--------------------------------------------------------
kpAddBaseNumbers(kp)
kpAddCytobandLabels(kp,force="all",cex=0.4)
#--------------------------------------------------------

#--------------------------------------------------------
# Data Panel 2 ###
# bam coverage plot data ###
# medium regions and their coverage
#--------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
#--------------------------------------------------------
gr.windows <- tileGenome(seqinfo(Hsapiens), tilewidth=100000,cut.last.tile.in.chrom=TRUE)
gr.data    <- GRanges(c("chr1", "chr2"), IRanges(c(10, 50), c(15, 55)), value = c(20, 10))
gr.data.RleList <-  mcolAsRleList(gr.data, varname = "value")
seqlevels(gr.windows, force=TRUE) <- names(gr.data.RleList)
gr.data.binnedAvg <- binnedAverage(gr.windows, gr.data.RleList, "value")
#--------------------------------------------------------
kpText(kp, chr=seqlevels(kp$genome), 
       y=0.4, x=0, data.panel = 2, r0=0.2, r1=0, 
       col="#444444", label="30x", cex=0.8, pos=2)
kpAbline(kp, 
         h=0.4,
         data.panel = 2, 
         r0=0.2,
         r1=0,
         col=dp.colors[3])
kpText(kp, chr=seqlevels(kp$genome), 
       y=0.4, x=0, data.panel = 2, r0=0.2, r1=0, 
       col="#444444", label="30x", cex=0.8, pos=2)
#--------------------------------------------------------
kpPlotCoverage(kp, data = mid.regs, r0=0.2, r1=0, col=dp.colors[2], data.panel = 2)
kpPlotCoverage(kp, data = mid.regs, r0=0.2, r1=0.12, col=dp.colors[1], data.panel = 2)
#--------------------------------------------------------
#
#--------------------------------------------------------
kpAbline(kp, 
         h=0.4,
         data.panel = 2, 
         r0=0.2,
         r1=0,
         col=dp.colors[3])
#--------------------------------------------------------

graphics.off()
#--------------------------------------------------------
q('no')




