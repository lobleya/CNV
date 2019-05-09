#-----------------------------------------------------------------
library("karyoploteR")
library("regioneR")
library("zoo")
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
#  karyoplote R: add in seqdata depth from normals
#                check the genomic binned ranges
#                check top 10 upreg and downreg genes
#--------------------------------------------------------
parser$add_argument('--logR', dest='logRfile', type="character",
                    help="logR data file for processing",
                    default="H://projects/ctDNA_pancreas/test.txt")
parser$add_argument('--tum', dest='Tfile', type="character",
                    help="logR data file for processing",
                    default="H://projects/ctDNA_pancreas/tum.txt")
parser$add_argument('--norm', dest='Nfile', type="character",
                    help="logR data file for processing",
                    default="H://projects/ctDNA_pancreas/norm.txt")
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
d<-d[which(d$depth>0),]
d$pvalue<-1-pnorm(abs(d$log2-mean(d$log2))/sd(d$log2))
d$fdr<-p.adjust(d$pvalue)
#--------------------------------------------------------
gmax<-2
chrs<-parser$chr
#--------------------------------------------------------
for(chr in paste("chr",1:22,sep=""))
{
#--------------------------------------------------------
dp.colors  <- c("#FFBD07", "#00A6ED",  "#FF3366", "#8EA604", "#C200FB")
#-----------------------------------------------------------------
data.ranges<- GRanges(IRanges(start=d$start,end=d$end),
                      seqnames=d$chromosome,
                      gene=d$gene,
			    gc=d$gc,
                      depth=d$depth,
                      logR=d$log2,
			    fdr=d$fdr,
			    y=d$log2,
                      pvalue=d$pvalue)
#--------------------------------------------------------
signif.points<-data.ranges[which(data.ranges$pvalue<0.05 & data.ranges@seqnames==chr & abs(data.ranges$logR)>1),]
smin<-round(max(abs(signif.points$logR))*-1)
smax<-round(max(signif.points$logR))
#-----------------------------------------------------------------
data.points.colors<-dp.colors
#-----------------------------------------------------------------
pdf(paste(chr,".pdf",sep=""),width=20,height=10)
#-----------------------------------------------------------------
gains <-d[which(d$log2>0),]
losses<-d[which(d$log2<0),]
#-----------------------------------------------------------------
gains.genes<-unique(gains$gene[which(gains$fdr<0.1)])
gains.genes.data<-data.ranges[match(gains.genes,d$gene)]
loss.genes<-unique(losses$gene[which(losses$fdr<0.1)])
loss.genes.data<-data.ranges[match(loss.genes,d$gene)]
#--------------------------------------------------------
mid.regs <- data.ranges[which(data.ranges$pvalue >= 0.05 & data.ranges@seqnames==chr),]
#--------------------------------------------------------
kp<-plotKaryotype(genome="hg38",
                  chromosomes=chr,
                  plot.type=2)
#-----------------------------------------------------------------
#   Data Panel 1                                                 #
#-----------------------------------------------------------------
#Big regions
#-----------------------------------------------------------------
big.regs.up<-data.ranges[which(data.ranges$pvalue <0.05 & data.ranges@seqnames==chr & data.ranges$log2>0 ),]
big.regs.up$x0<-big.regs.up@ranges@start
big.regs.up$x1<-big.regs.up$x0+big.regs.up@ranges@width
#-----------------------------------------------------------------
big.regs.down<-data.ranges[which(data.ranges$pvalue <0.05 & data.ranges@seqnames==chr & data.ranges$log2 < 0),]
big.regs.down$x0<-big.regs.down@ranges@start
big.regs.down$x1<-big.regs.down$x0+big.regs.down@ranges@width
#--------------------------------------------------------
kpRect(kp, data = mid.regs, 
        col="lightgray", 
        y0=0.4,
        y1=0.8,
       border="lightgray", 
        r0=abs((smax/smin)), 
       data.panel = 1)
#--------------------------------------------------------
#Data points
yrange<-max(data.ranges$y)-min(data.ranges$y)
#--------------------------------------------------------
ymin<-  5
ymax<- -5
N<-10
#--------------------------------------------------------
kpAxis(kp, 
       ymin = round(smin), 
       ymax = round(smax), 
       r0=0.2, 
       r1=1, 
       numticks = smax-smin+1, 
       col="#666666", cex=0.5)
#--------------------------------------------------------
yh<-1/yrange*abs(min(data.ranges$y))
#--------------------------------------------------------
#Mean and sd of the data points.  
chr.dp <- sort(keepSeqlevels(x = data.ranges, 
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
        col =dp.colors[3])
#--------------------------------------------------------
kpPlotRibbon(kp, 
             chr=chr, 
             data=chr.dp[3:(length(chr.dp)-3)],
             y0=rmean-rsd,
             y1=rmean+rsd,
             ymin=min(rmean-rsd),
             ymax=max(rmean+rsd),
             r0=abs(smax/smin), 
             r1=1,
             col="#FF336633", 
             border=NA)
#--------------------------------------------------------
kpPoints(kp, 
         data=signif.points,
         y=signif.points$log2,
         pch=16, 
         cex=0.5, 
         col="brown", 
         r0=0.2, 
         r1=1, 
         ymin=smin,
         ymax=smax)
#--------------------------------------------------------
ymax=max(rmean)
ymin=min(rmean)
#--------------------------------------------------------
kpText(kp,chr=chr,y=0.4,x=0,data.panel=1,ymin=0.2,ymax=0.4,label="logR",
       r0=1,r1=1.2,pos=2,cex=0.8,col="black")

#--------------------------------------------------------
signif.points<-signif.points[order(signif.points$pvalue),]
signif.genes<-signif.points[which(signif.points$gene!="-"),]
UN<-which(signif.genes$y>0)
DN<-which(signif.genes$y<0)
up.signif.genes<-signif.genes[UN[1:min(N,length(UN))],]
down.signif.genes<-signif.genes[DN[1:min(N,length(DN))],]
#----------------------------------------------------------
top.genes    <-reduce(split(signif.genes, 
                            elementMetadata(signif.genes)$gene),
                            drop.empty.ranges=TRUE)
top.genes    <-top.genes[which(lapply(top.genes,length)>0)]
                      
top.pvals    <-lapply(split(signif.genes,
                               signif.genes$gene),function(x){return(min(x$pvalue,na.rm=T))})

top.pvals    <-order(unlist(top.pvals[match(names(top.genes),
                                            names(top.pvals))]))
#---------------------------------------------------------
gene.names<-names(top.genes)[top.pvals]
top.genes <-top.genes[top.pvals]
#---------------------------------------------------------
ug<-which(names(top.genes) %in% up.signif.genes$gene)
dg<-which(names(top.genes) %in% down.signif.genes$gene)
#---------------------------------------------------------
up.genes  <-top.genes[ug, ]
down.genes<-top.genes[dg, ]
#--------------------------------------------------------
if(length(up.genes)>1)
{
kpPlotMarkers(kp,
              chr=chr,
              x=unlist(lapply(up.genes,function(x){return(min(ranges(x)@start))})),
              labels = gsub(",.+$","",names(up.genes)), 
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
              labels = gsub(",.+$","",names(down.genes)), 
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
#library(BSgenome.Hsapiens.UCSC.hg38)
#library(data.table)
#--------------------------------------------------------
#
#--------------------------------------------------------
gr.windows        <-  tileGenome(seqinfo(Hsapiens)[chr], 
                                 tilewidth=50000,
                                 cut.last.tile.in.chrom=TRUE)
gr.data           <-  data.ranges[which(data.ranges@seqnames==chr),]
#---------------------------------------------------------
kpPlotDensity(kp,gr.data,data.panel=2,col="yellow",border="gold",r1=0.3,ymin=700,ymax=0)
kpPlotDensity(kp,gr.data,data.panel=2,border="royalblue",col="lightblue",r1=0.3,ymin=1000,ymax=0)
#---------------------------------------------------------
kpText(kp,chr=chr,y=0.4,x=0,data.panel=2,ymin=1000,ymax=0,label="30x ",
       r0=0.2,r1=0.3,pos=2,cex=0.8,col="black")
kpAbline(kp,h=30,data.panel=2,ymin=700,ymax=0,r1=0.3,col="black")

graphics.off()
}
#--------------------------------------------------------



