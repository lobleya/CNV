library("karyoploteR")
#--------------------------------------------------------
d<-read.table("C://Users/lobley01/Documents/ctDNA_PROJECT/test.txt",
              header=T)
gmax<-2
chrs<-"chr9"

gains<-d[which(d$log2>0),]
losses<-d[which(d$log2<0),]


data.ranges<-GRanges(chr=d$chr,
                     IRanges(start=d$start,end=d$end),
                     seqnames=d$chr,genes=d$gene,
                     logR=d$log2,pvalue=d$pvalue)

gains.genes<-unique(gains$gene[which(gains$pvalue<0.05)])
gains.genes.data<-data.ranges[match(gains.genes,d$gene)]
loss.genes<-unique(losses$gene[which(losses$pvalue<0.05)])
loss.genes.data<-data.ranges[match(loss.genes,d$gene)]
data.ranges<-toGranges(d)
kp<-plotKaryotype(genome="hg38",chromosomes=c(chrs),plot.type=1)
#--------------------------------------------------------
data.ranges$depth<-d$depth
kp.regs<-data.ranges
y<-d$log2
y<-ifelse(d$log2>1.5,1.5,ifelse(d$log2< -1.5,-1.5,d$log2))
#--------------------------------------------------------
kpHeatmap(kp, 
          data.ranges[seqnames(data.ranges)==chrs], 
          y=y, 
          ymax=1.5,
          ymin=-1.5,
          colors = c("blue", "white", "red"), 
          r0=0,r1=0.1)
#--------------------------------------------------------
kpPlotRegions(kp, 
              kp.regs, 
              r0 = 0.15, 
              r1 = 0.15+0.001*data.ranges$depth, 
              col="navy")
#--------------------------------------------------------
kpAddCytobandLabels(kp,
                    cex=0.5)
#--------------------------------------------------------
kpAddBaseNumbers(kp, 
                 tick.dist=10000000, 
                 minor.tick.dist=1000000)
#--------------------------------------------------------
kpAddCytobandLabels(kp)
#--------------------------------------------------------
#show signif genes only i.e. zoom in on regions
#--------------------------------------------------------


