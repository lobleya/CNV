#---------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(karyoploteR)
#---------------------------------------------------------------
blues<-colorRampPalette(c("white","royalblue"))(10)
reds<-colorRampPalette(c("white","brown"))(10)
greens<-colorRampPalette(c("white","forestgreen"))(10)
#---------------------------------------------------------------
cyto.cols      <-getCytobandColors()
cyto.bands     <-getCytobands()
cyto.bands$cols<-cyto.cols[match(as.character(cyto.bands$gieStain),as.character(names(cyto.cols)))]
#---------------------------------------------------------------
# Add colored rectangles for chrome segments
# Add in labels for sigdiff segments with p-values
#---------------------------------------------------------------
#---assign bands to cnr ----------------------------------------
cyto.bands<-as.data.frame(cyto.bands)
CNS<-lapply(CNS,function(x){x$karyo.id<-0; return(x)})
#---------------------------------------------------------------
col2rgb<-function(hex,alpha){ rcol<-col2rgb(hex,alpha=T)[,1]; 
                              return(rgb(rcol[1],rcol[2],rcol[3],alpha=alpha,maxColorValue = 255))     
}
#---------------------------------------------------------------
cleanCNR<-function(CNR)
{
  CNR<-lapply(CNR,function(cnr){
    cnr<-cnr[which(cnr$depth > 1),]
    cnr<-cnr[which(cnr$gene!="Antitarget"),]
    cnr<-cnr[grep("chr[XY]",cnr$chromosome,invert=T),]
    return(cnr)
  })
  
  return(CNR)
}
#---------------------------------------------------------------
cleanCNS<-function(CNS)
{
  CNS<-lapply(CNS,function(cns){
    tmp=unlist(strsplit(as.character(cns$gene),split=","))
    #---------------------------------------------
    tmp=tmp[grep("ENSE",tmp,invert=T)]
    tmp=tmp[grep("\\-[23]",tmp,invert=T)]
    cns$gene=paste(unique(tmp),collapse=",")
    #---------------------------------------------------
    return(cns)
  })  
  return(CNS)
}
#---------------------------------------------------------------
# clean CNS and CNR
#---------------------------------------------------------------
CNS<-cleanCNS(CNS)
CNR<-cleanCNR(CNR)
#---------------------------------------------------------------
#---- for karyo plots ----------------

#---- for chr plots ------------------
setCytoBands<-function(CNS,cyto.bands)
{
 CNS<-lapply(CNS,function(cns){
  
  for(i in 1:nrow(cyto.bands))
  {
    w=which(as.character(cns$chromosome)==as.character(cyto.bands$seqnames[i]) & 
              cns$start >= cyto.bands$start[i] & cns$start <= cyto.bands$start[i] + cyto.bands$width[i] )
    cns$karyo.id[w]=i
  }
  return(cns)
 })
}
#---------------------------------------------------------------
setSegMeans<-function(CNR,CNS,ARM)
{
   #-------------------------
   j=1
   #-------------------------
   for(j in 1:length(CNR))
   {
    cns<-CNS[[j]]
    cnr<-CNR[[j]]
    cnr$cns.id=0
    cnr$karyoband.id=0
    cnr$cols="white"
    cnr$cns.mean=0
    cnr$cns.sd=0
    cnr$astart=1:nrow(cnr)
    
 
    for(i in 1:nrow(cns))
    {
     w<-which(cnr$start >= as.numeric(cns$start[i]) & cnr$start <= as.numeric(cns$end[i]) & 
			as.character(cnr$chromosome)==as.character(cns$chromosome[i]))
     cnr$cns.id[w]  <-i  
     cnr$cns.mean[w]<-weighted.mean(cnr$log2[w],weights=cnr$weight[w],trim=0.95,na.rm=T)
     cnr$cns.sd[w]  <-(cnr$log2[w]-cnr$cns.mean[w][1])/sd(cnr$log2[w])
    }
 
    for(i in 1:nrow(cyto.bands))
    {
     w<-which(cnr$start >= as.numeric(cyto.bands$start[i]) & 
              cnr$start <= as.numeric(cyto.bands$start[i]+cyto.bands$width[i]) & 
                as.character(cnr$chromosome)==as.character(cyto.bands$seqnames[i]))
              
     cnr$karyoband.id[w]  <-i  
     cnr$cols[w]          <-cyto.bands$col[i]
    }

    
    cnr$karyo.mean=0
    cnr$karyo.sd=0
    cnr$karyo.label=0
    cnr$karyo.id=0
    cnr$karyo.col=0
    
    for(i in 1:nrow(ARM))
    {
      w<-which(cnr$start >= ARM$start[i] & 
                 cnr$start <= ARM$end[i] & 
                  as.character(cnr$chromosome)==as.character(ARM$chr[i]))
      
      cnr$karyo.mean[w]<-weighted.mean(cnr$log2[w],cnr$weight[w],trim=0.95)
      cnr$karyo.sd[w]  <-cnr$log2[w]-cnr$karyo.mean[w]
      cnr$karyo.id[w]  <-i  
      cnr$karyo.col[w] <-ARM$col[i]
      
    }
    
    CNR[[j]]=cnr
  
   }
 
  
}
#----------------------- plot by chromosomes --------------------
#
# w=which(cnr$chromosome %not in% c("chrX","chrY"))
#----------------------------------------------------------------
storelposx1=c()
storelposx2=c()
#----------------------------------------------------------------
#
#----------------------------------------------------------------
makeChrPlot<-function(cnr,cns,chr,cname)
{
  #-------------------------------------------------
  par(mar=c(3,2,3,1),las=2)
  w=which(cnr$chromosome==chr & cnr$cns.sd!=0)
  #-------------------------------------------------
  plot(cnr$astart[w][which(abs(cnr$cns.sd[w])<0.2)], 
       cnr$log2[w][which(abs(cnr$cns.sd[w])<0.2)],
       type="p",
       xaxt="n",
       pch =20,
       ylim=c(-1.5,1.5),
       xlab="",
       ylab="logR",
       col =blues[cnr$weight[w]*10])
  #-------------------------------------------------
  title(paste(cname,sep=" "),adj=0,line=0.5)
  #--------------------------------------------------
  # Shading add chr band colors to plots
  #--------------------------------------------------
  for(k in unique(cnr$karyoband.id[w]))
  {
    rect(min(cnr$astart[w][which(cnr$karyoband.id[w]==k)]),
         -1.5,
         max(cnr$astart[w][which(cnr$karyoband.id[w]==k)]),
         1.5,
         col=cyto.bands$cols[k],
         border="gray80",
         density=40)
  }
  #--------------------------------------------------
  # Add points to plots i.e. mean points
  #--------------------------------------------------
  points(cnr$astart[w][which(abs(cnr$cns.sd[w])<0.5)],
         cnr$log2[w][which(abs(cnr$cns.sd[w])<0.5)],
         col=blues[cnr$weight[w]*60],pch=20)
  points(cnr$astart[w],cnr$cns.mean[w],ylim=c(-2,2),col="orangered",pch=20)
  #--------------------------------------------------
  # Add labels to plots i.e. signif logR segs only
  #--------------------------------------------------
  smeans<-cnr$cns.mean/sd(cnr$cns.mean)
  sidx  <-unique(cnr$cns.id[w][which(1-pnorm(abs(smeans[w]))<0.05)])
  sup   <-which(cns$cns.mean[w][sidx]>0)
  sdown <-which(cns$cns.mean[w][sidx]<0)
  if(length(sup)>5)   { o=order(cns$cns.mean[w][sup],decreasing=T); sup=sup[o[1:5]]  }
  if(length(sdown)>5) { o=order(cns$cns.mean[w][sdown],increasing=T); sdown=sdown[o[1:5]]  }
  sidx<-c(sup,sdown)
  #--------------------------------------------------
  for(k in sidx)
  {
   mini<- min(cnr$astart[which(cnr$cns.id==k)])
   maxi<- max(cnr$astart[which(cnr$cns.id==k)])
   top <- ifelse(cns$log2[k]>0, 3, -1.8) #assign label to top or bottom
   xx  <- cnr$astart[mini]
   #--- alternate labels above and below plots
   lgene<-cns$gene[k]
   lmax <-strwidth(cyto.bands$name[k])
   #----------------------------------------------
   if(length(lgene)>0)
   { 
     lgene<-unlist(strsplit(as.character(lgene),","));
     lgene<-lgene[grep("^AP0",lgene,invert=T)]
     lgene<-lgene[1:min(15,length(lgene))]
     lmax<-max(unlist(sapply(lgene,strwidth)))
   }
   #----------------------------------------------
   # draw lgene labels
   #----------------------------------------------
   if(which(k==sidx)==1)
   {
    lposx1<-xx-100
    lposx2<-xx+lmax*0.8
   }else{
     lposx1<-lposx1+200
     lposx2<-lposx1+lmax*0.8
   }
   #----------------------------------------------
   storelposx1<-c(storelposx1,lposx1)
   storelposx2<-c(storelposx2,lposx2)
   #----------------------------------------------
   rect(xleft=lposx1,top-0.1,
        xright=lposx2,top-0.23*max(1,length(lgene))-0.2,
        col="beige",border="orange",lwd=1,xpd=T)
   #----------------------------------------------
   top=top-0.2
   text(x=lposx1+100,
        y=top,
        xpd=T,
        bg="beige",
        cex=0.8,
        adj=0,
        cyto.bands$name[cns$karyo.id[k]])
   #------------------------------------------------
   if(length(lgene)>0)
   {
     top=top-0.05
     for(l in 1:length(lgene))
     {
      top=top-0.1
        
      text(x=xx,
           y=top,
           cex=0.6,
           xpd=T,
           adj=0,
           lgene[l])
       
     }
   }
   
   #--- finish results ---------------------------
  }#for all signif regions  
  #-------------------------------------------------
  abline(h=-1,col="pink")
  abline(h=0, col="gray") 
  abline(h=1, col="pink")
#-----------------------------------------------------------------
}
#-------------------------------------------------
plotAllChr<-function(cns,cnr)
{
  
  #-----------------------------------------#
  # plotAllChr
  #-----------------------------------------#
  w=which(abs(cnr$karyo.sd)<0.01 & cnr$karyo.mean!=0)
  #---- plot mean by KARYO P|Q SEGMENTS ----#
  zscores=cnr$karyo.mean/sd(cnr$karyo.mean)
  signif=1-pnorm(abs(zscores))
  sidx  =unique(cnr$karyo.id[which(signif<0.1)])
  #------------------------------------------
  # signif colors
  #------------------------------------------
  cnr$karyo.col=blues[1+as.integer(cnr$weight*30)]
  up=which(signif<0.05 & zscores >0)
  down=which(signif<0.05 & zscores <0)
  cnr$karyo.col[up]=reds[1+as.integer(cnr$weight[up]*30)]
  cnr$karyo.col[down]=greens[1+as.integer(cnr$weight[down]*30)]
  #------------------------------------------
  cnr$mcolor="black"
  cnr$mcolor[intersect(signif<0.05,zscores>0)]="red"
  cnr$mcolor[intersect(signif<0.05,zscores<0)]="forestgreen"
  sx   = seq(1,44,by=2)
  sx1  = unlist(sapply(sx,function(x){return(which(cnr$karyo.id==x)[1])}))
  sx2  = unlist(sapply(sx,function(x){return(max(which(cnr$karyo.id==x)))}))
  shade= data.frame(sx1=sx1,sx2=sx2,sy1=rep(-0.5,length(sx1)),sy2=rep(0.5,length(sx2)))
  xpos = 200+cnr$astart[sx1]
  glab = data.frame(gt=unique(ARM[,2]),xpos=xpos,ypos=rep(c(-0.25,-0.2),11))
  #-----------------------------------------
  g=ggplot(cnr[which(cnr$cns.sd!=0),],aes(x=astart,y=karyo.mean))+ # should be 300 different values here  
    geom_point(data=cnr[w,],
               aes(x=cnr$astart[w],y=cnr$log2[w]),
               col=cnr$karyo.col[w],size=0.1)+   # smoothed points
        geom_point(data=cnr,aes(x=cnr$astart,y=cnr$karyo.mean),col=cnr$mcolor,size=0.5)+     # smoothed means for karyo bins
            geom_hline(yintercept=0.1,col="pink")+ 
             geom_hline(yintercept=0,col="#cccccc")+
              geom_hline(yintercept=-0.1,col="pink")+
                geom_vline(xintercept=cnr$astart[sx1],size=0.2, col="navy", linetype="dashed",alpha=0.5)+
                 geom_vline(xintercept=cnr$astart[sx2],size=0.2, col="navy", linetype="dashed",alpha=0.5)+ 
                  geom_label(inherit.aes=F,size=3,data=glab,aes(x=xpos, y=ypos, label=gt))+
			ylim(c(-0.3,0.3))+ylab("logR")+xlab("genome") # add logRs
  
  #------------------------------------------#
  return(g)
  #------------------------------------------#
}
#-----------------------------------------------------------------
#end of crappy functions !!!
#-----------------------------------------------------------------