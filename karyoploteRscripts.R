library("karyoploteR")
#--------------------------------------------------------
gmax<-2
chrs<-"chr9"
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
#show losses and gains

#show signif genes only i.e. zoom in on regions


#data <- toGRanges(data.frame(chr="chr1", start=10000000*(0:23), end=10000000*(1:24)))

#kp <- plotKaryotype("hg19", plot.type=2, chromosomes=c("chr1", "chr2"))

#We can specify all data values separately. If missing y0, it defaults to ymin
#kpBars(kp, chr=as.character(seqnames(data)), x0=start(data), x1=end(data), y1=y1,
#       col="#FFBBBB", border="#EEAAAA")
#kpLines(kp, data=data, y=y1, col="red")

#or we can provide all data into a single GRanges object
#mcols(data) <- data.frame(y0=y0, y1=y1)
#kpBars(kp, data[data$y0>data$y1], col="orange", border="orange", data.panel=2)
#kpBars(kp, data[data$y0<=data$y1], col="purple", border="purple", data.panel=2)

#kpLines(kp, data, y=data$y1, data.panel=2, col="red")
#kpLines(kp, data, y=data$y0, data.panel=2, col="blue")

#kpAxis(kp, data.panel = 1, cex=0.8, numticks = 5, col="#777777")
#kpAxis(kp, data.panel = 2, cex=0.8, numticks = 5, col="#777777")

