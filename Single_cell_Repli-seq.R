
####Author: Vishnu Dileep
####Updated: Jan 2018

###see library 'copynumber' for parameter descriptions
###Function to Segment genome based on copynumber####
segment_vd=function(data,my_gamma=3,my_kmin=10,my_fast=F){

  pcf.seg=NULL
  for(chr in 1:20){
    cat("\nCHR=",chr)
    segmentation=NULL
    datab=subset(data,CHR==chr)
    datab[,2:ncol(datab)]=na.approx(datab[,2:ncol(datab)],na.rm=F)
    datab[,2:ncol(datab)]=na.locf(datab[,2:ncol(datab)],na.rm=F)
    datab[,2:ncol(datab)]=na.locf(datab[,2:ncol(datab)],na.rm=F,fromLast = T)
    datab$CHR=NULL
    
    segmentation=pcfPlain(datab,gamma=my_gamma,kmin=my_kmin,return.est = F,fast=my_fast)
    segmentation$chrom=chr
    pcf.seg=rbind(pcf.seg,segmentation)
  }
  pcf.seg=pcf.seg[order(pcf.seg$sampleID,pcf.seg$chrom,pcf.seg$start.pos),]
  return(pcf.seg)
}
###########################################################################
##########################Load libraries#####################################
###########################################################################

source("./scRT_functions.R")

library(copynumber)
library(zoo)
library(circlize)
library(gridExtra)
library(lattice)
library(ggplot2)
library(mixtools)
library(fBasics)


###########################################################################
##########################Read datasets#####################################
###########################################################################

#######Read and Pre-process single cell data (haploid cells)###

RT=read.delim("Hap_96cells_mm10_named.txt", header=T)   ###

read=stack(colSums(RT[,3:ncol(RT)]));names(read)=c("mapped_reads","cell")
remove=read[read$mapped_reads<=250000,]  ##List of cells with below 250,000 reads
remove=as.character(remove$cell)
remove=c(remove,"C5","C8","G1","G6","H2","H4","H11","H12") ##Add list of cell with ploidy abberations (identifed visually by plotting all chromosomes)##
RT=RT[,!names(RT) %in% remove]

for(i in 3:ncol(RT)) RT[,i]=RT[,i]* (1000000/sum(RT[,i])) ###Convert to Reads per Million
RT=RT[RT$CHR!="chrM" & RT$CHR!="chrY",]  ##Remove Mitochindrila DNA and Y chrom##

RT[,1]=as.character(RT[,1])   ##Re-format data for copynumber package###
RT[RT$CHR=="chrX",1]="chr20"
RT$CHR=as.numeric(substring(RT$CHR,4,6))

##########read Population RT dataset (bulk)###
RTbulk=read.delim("./cas.129.bulkRT.txt", header=T)
names(RTbulk)[3]="cast.129"
RTbulk=RTbulk[RTbulk$CHR!="chrM" & RTbulk$CHR!="chrY",]

RTbulk[,1]=as.character(RTbulk[,1])     ##Re-format data for copynumber package###
RTbulk[RTbulk$CHR=="chrX",1]="chr20"
RTbulk$CHR=as.numeric(substring(RTbulk$CHR,4,6))


#####Aliging Bulk RT and single-cell RT##
RT=merge(RT,RTbulk,by=c("CHR","POS"))
RT=RT[order(RT$CHR,RT$POS),]
RTc=RT
RTc$cast.129=NULL

###########################################################################
###Identify and remove bad bins and bad segments based on G1 and G2 cells##
###########################################################################  
g_con=RTc[,c("CHR","POS","H5","H6","H8","H9")]  ##These cells were sorted for G1/G2 ####
g_con$mean=rowMeans(g_con[,3:6])
g_con=g_con[,c("CHR","POS","mean")]



##find bad 50Kb bins based on G1/G2 data###
badwins=subset(g_con, mean<=quantile(g_con$mean,0.05) | mean>=quantile(g_con$mean,0.99) )
badwins=badwins[,c("CHR","POS")]
badwins$ID=as.character(paste(badwins$CHR,"_",badwins$POS,sep=""))

###Find bad segments based on G1/G2 data###
gseg=segment_vd(g_con[,c("CHR","POS","mean")])
gsegframe=interpolate.pcf(gseg,g_con[,1:2])
names(gsegframe)=c("CHR","POS","mean")

goodsegs=subset(gsegframe,mean>quantile(gseg$mean,0.05))
goodsegs$mean=NULL

###Mask outlier bins regions with NA in single cell dataset####
RTc$ID2=as.character(paste(RTc$CHR,"_",RTc$POS,sep=""))
L=RTc$ID2 %in% badwins$ID
RTc[L,3:ncol(RTc)]=NA
RTc$ID2=NULL
###Remove bad segments from single cell dataset###
RTc=merge(RTc,goodsegs,by=c("CHR","POS"))
RTc=RTc[order(RTc$CHR,RTc$POS),]

###########################################################################
###Normalize all single-cell data with average G1/G2 single cells##########
###########################################################################
  base=(RTc$H5 + RTc$H6 + RTc$H8 + RTc$H9)/4
  RTc$base=base

for(i in 3:ncol(RTc)) RTc[,i]=RTc[,i]/RTc$base
RTc$base=NULL 

RTc[,3:ncol(RTc)]=do.call(data.frame,lapply(RTc[,3:ncol(RTc)], function(x) replace(x, !is.finite(x),NA))) #replace nonfinite with NA


#####Center and scale data#####

for(i in 3:ncol(RTc)) RTc[,i]= scale(RTc[,i],scale=FALSE)
RTc=rescaleDatasets(RTc)

###########################################################################
########################remove control cells###########################
###########################################################################

remove=c("H5","H6","H8","H9")
RTc=RTc[,!names(RTc) %in% remove]

########Perform median filtering one CHR at a time#######
chrs=levels(as.factor(as.character(RTc$CHR)))
AllData=NULL
for(chr in chrs){
cat("\n",chr,"\n")
chrData=RTc[RTc$CHR==chr,1:2]
full=RTc[RTc$CHR==chr,1:2]
full$ID2=as.character(paste(full$CHR,"_",full$POS,sep=""))
for(i in 3:ncol(RTc))
{
cat(" ",names(RTc)[i]," ")
med.temp=RTc[RTc$CHR==chr,c(1,2,i)]
med.temp=na.omit(med.temp) 
up.lim=quantile(med.temp[,3],0.99,na.rm=T)
down.lim=quantile(med.temp[,3],0.01,na.rm=T)
med.temp[med.temp[,3]>=up.lim | med.temp[,3]<=down.lim,3]=NA
med.temp=na.omit(med.temp) 
med.temp$smooth=runmed(med.temp[,3],15) 
med.temp$ID=as.character(paste(med.temp$CHR,"_",med.temp$POS,sep=""))
L=!full$ID2 %in% med.temp$ID
med.temp$ID=NULL
omitted=full[L,c(1,2)]
add=cbind(omitted[,1:2], matrix(data = NA, nrow = nrow(omitted), ncol = 2))
names(add)=names(med.temp)
med.temp=rbind(med.temp,add)
med.temp=med.temp[order(med.temp$CHR,med.temp$POS),]
chrData=cbind(chrData,med.temp[,4])
}
AllData=rbind(AllData,chrData)
}
names(AllData)=names(RTc)
 
RTs=AllData

##########################################################################
##########################Segmentation using copynumber package #####################################
###########################################################################

###Seg each chromosome using pcfPlain##
pcf.seg=segment_vd(RTs,my_fast = F,my_gamma = 3,my_kmin = 5)
segframe=interpolate.pcf(pcf.seg,RTs[,1:2])
names(segframe)[1:2]=c("CHR","POS")

#check.seg(RTs,segframe,sample=c("B9","B3","B2","A11"),binary=F,chrs=c(1,2,16),x1=2e6,x2=200e6) ###check segmentation by plotting##

###########################################################################
##########################Binarize data#####################################
###########################################################################

copyframe=segframe

plotmix=function(linep){     ####function to plot results of mixture model###
  d=density(test,adjust=0.5)
  xa=min(d$x)-0.1; xb=max(d$x)+0.1
  ya=0; yb=max(d$y)+0.1
  
  plot(d,xlim=c(xa,xb),ylim=c(ya,yb))
  abline(v=linep)
  par(new=T)
  d=density(rnorm(1000000, mean=model$mu[1],sd=model$sigma[1]),adjust=0.5)
  plot(d$y*model$lambda[1]~d$x,type="l",col="red",xlim=c(xa,xb),ylim=c(ya,yb))
  d=density(rnorm(1000000, mean=model$mu[2],sd=model$sigma[2]),adjust=0.5)
  lines(d$y*model$lambda[2]~d$x,col="green",type="l")
  legend("topright",paste(sample,round(mean.diff,2),round(skew,2),sep=" "),bty="n")
}

#####Running 100 different thresholds to find best threshold for binarization for each single-cell####

testing=NULL
c.noskew=NULL;c.lskew=NULL;c.rskew=NULL;
samples=names(copyframe)[3:ncol(copyframe)]       ##list of single-cells
point2Cs=NULL                                     ###vector to store the best threshold for each cell##
for(sample in samples){                           ##Iterate over each sample###
  seg=NULL;model=NULL
  cat("\n",sample,"\n")
  test=segframe[,which(names(segframe)==sample)]
  test=test[abs(test)<2]                          ###mask bins with extreme segmented values for binarization threshold calculation###
  model <- normalmixEM(x=test, k=2)               ###Mixture model with 2 components (Replicated and un-replicated)
  
  skew=skewness(test)                            ###Mixture models fails in cells in very early or late is S-pahse. Calculating
  mean.diff=abs(diff(model$mu))                  ###skew in the data and difference of component means in the mixture model to 
  cat("\n","skew - ",skew)                       ###identify these cells.
  cat("\n","Diff - ",mean.diff)
  
  testing=rbind(testing,data.frame(sample,mean.diff,skew))
  ###Setting range of the thresholds to iterate over##
  sweep.min=min(test)
  sweep.max=max(test)
  
  cpoints=seq(sweep.min,sweep.max,(sweep.max-sweep.min)/100)  ##Thresholds to iterate over###
  cpoints.cors=NULL
  seg=data.frame(test=test,testb=NA,testc=NA)
  
  if(mean.diff>0.7){       ###If difference in component means are above 0.7 (may need to be adjusted emperically) then there is
                           ##a clear replicated and un-replicated fraction in the data####
    plotmix(model$mu)
    for(cpoint in cpoints ){     ##Iterate over all Thresholds
      seg[,"testb"]=seg[,"test"]-cpoint        ###subtract threshold from data##
      seg[seg$testb<0,"testc"]=min(model$mu)   ###Binarize with component means instead of 1s and 0s
      seg[seg$testb>0,"testc"]=max(model$mu)
      cpoints.cors=c(cpoints.cors,sum(abs(seg$test-seg$testc)))  ##Calculate Mahattan distance between Binarized and Un-brinarized data
    }
                     }
  
  
  
  if(mean.diff<0.7){     ###If difference in compmonent means are above 0.7 (maybe need to be adjusted emperically) then the cells is most
                         ###likely in very early or very late S-phase. Detemrine based on skew of data##
    if(skew<(-0.2)){
      c.lskew=c(c.lskew,sample)
      cat("\n","skew to left-late S")
      plotmix(quantile(test,c(0.05,0.5)))
      for(cpoint in cpoints ){
        seg[,"testb"]=seg[,"test"]-cpoint
        seg[seg$testb<0,"testc"]=quantile(test,0.05)
        seg[seg$testb>0,"testc"]=quantile(test,0.5)
        cpoints.cors=c(cpoints.cors,sum(abs(seg$test-seg$testc)))
      }
    }
    
    if(skew>0.2){
      c.rskew=c(c.rskew,sample)
      cat("\n","skew to right-early S")
      plotmix(quantile(test,c(0.5,0.95)))
      for(cpoint in cpoints ){
        seg[,"testb"]=seg[,"test"]-cpoint
        seg[seg$testb<0,"testc"]=quantile(test,0.5)
        seg[seg$testb>0,"testc"]=quantile(test,0.95)
        cpoints.cors=c(cpoints.cors,sum(abs(seg$test-seg$testc)))
      }
    }
    
    
    if(skew>=(-0.2) & skew<=0.2) {          ###These cell may have noisy data, consider discarding later ###
      c.noskew=c(c.noskew,sample)
      cat("\n","no skew")
      plotmix(quantile(test,c(0.25,0.75)))
      for(cpoint in cpoints ){
        seg[,"testb"]=seg[,"test"]-cpoint
        seg[seg$testb<0,"testc"]=quantile(test,0.25)
        seg[seg$testb>0,"testc"]=quantile(test,0.75)
        cpoints.cors=c(cpoints.cors,sum(abs(seg$test-seg$testc)))
      }
    }
    
  }
  
  
  point.2C=cpoints[which(cpoints.cors==min(cpoints.cors,na.rm=T))]   ###Find the threshold at which the distance between binarized and un-binarized  is minimum
  point2Cs=c(point2Cs,point.2C[1])
# l=readline();if(l=="q") break
}
testing=testing[order(testing$mean.diff),]

###Binarize each single-cell using the specific ideal threshold calculated from above step and labels with 1s and 0s###

for(i in 3:ncol(copyframe)){
  cat("\n",names(copyframe)[i]," ",point2Cs[i-2],"\n")
  copyframe[,i]=copyframe[,i]-point2Cs[i-2]
  copyframe[copyframe[,i]<0,i]=0
  copyframe[copyframe[,i]>0,i]=1
}

#check.seg(RTs,copyframe,samples=c("B3","D2"),binary=T,chr=c(1,16),x1=2e6,x2=160e6)   ####check binarization, make binary=T in parameters###
###########################################################################
##########################Rank in S-phase#####################################
###########################################################################

rank=data.frame(samples=names(copyframe[,3:ncol(copyframe)]))
rankfun=function(x) length(x[x==1 & is.finite(x)])/length(x[is.finite(x)])*100
rank$rank=apply(copyframe[,3:ncol(copyframe)],2,rankfun)
rank=rank[order(rank$rank),]


#############FACS accuray calculation based on known sorting gates#########

rank2=rank[rank$samples %in% c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12"),]   ##early-S gate###
d=density(rank2$rank,adjust=0.5)
hist(rank2$rank,breaks=seq(0,100,5),freq=T,xlab="Percentage replicated",ylab="Cells",yaxt="n",xaxt="n",border="red",
     xlim=c(0,100),ylim=c(0,12), density=20,col="red")
axis(2, at=seq(0,30,2))
axis(1, at=seq(0,100,10))
par(new=T)
plot(d$y~d$x,type="l",col="red",xlim=c(0,100),ylim=c(0,0.3),yaxt="n",xaxt="n",ylab="",xlab="",lwd=2)



rank2=rank[!rank$samples %in% c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12","G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12"),]  ###Mid-gate, remove all early and late-S gate cells
d=density(rank2$rank,adjust=0.5)
par(new=T)
hist(rank2$rank,breaks=seq(0,100,5),freq=T,xlab="Percentage replicated",ylab="Cells",yaxt="n",xaxt="n",border="black",
     xlim=c(0,100),ylim=c(0,12), density=10,col="black")
axis(2, at=seq(0,30,2))
axis(1, at=seq(0,100,10))
par(new=T)
plot(d$y/10~d$x,type="l",col="black",xlim=c(0,100),ylim=c(0,0.006),yaxt="n",xaxt="n",ylab="",xlab="",lwd=2)



rank2=rank[rank$samples %in% c("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12"),] ##Late-S gate###
d=density(rank2$rank,adjust=0.5)
par(new=T)
hist(rank2$rank,breaks=seq(0,100,5),freq=T,xlab="Percentage replicated",ylab="Cells",yaxt="n",xaxt="n",border="blue",
     xlim=c(0,100),ylim=c(0,12), density=20,col="blue")
axis(2, at=seq(0,30,2))
axis(1, at=seq(0,100,10))
par(new=T)
plot(d$y~d$x,type="l",col="blue",xlim=c(0,100),ylim=c(0,0.3),yaxt="n",xaxt="n",ylab="",xlab="",lwd=2)


######################################################################
###################MErge RTc, copyframe to Bulk#################################
######################################################################
###merge bulk RT data and Binarized data##

copyframe=merge(copyframe,RTbulk,by=c("CHR","POS"))
copyframe=copyframe[order(copyframe$CHR,copyframe$POS),]

###merge bulk RT data and normalized scaled data##
RTs=merge(RTs,RTbulk,by=c("CHR","POS"))
RTs=RTs[order(RTs$CHR,RTs$POS),]

######################################################################
###################Add in all 50Kb bin POS#################################
######################################################################

####adding missing bins to binarized data###
copyframe$ID=as.character(paste(copyframe$CHR,"_",copyframe$POS,sep=""))


RT=read.delim("Hap_96cells_mm10_named.txt", header=T)
full=RT[,1:2]
names(full)=c("CHR","POS")
full=full[full$CHR!="chrM" & full$CHR!="chrY",]

full[,1]=as.character(full[,1])
full[full$CHR=="chrX",1]="chr20"
full$CHR=as.numeric(substring(full$CHR,4,6))
full$ID3=as.character(paste(full$CHR,"_",full$POS,sep=""))

L=!full$ID3 %in% copyframe$ID
copyframe$ID=NULL

omitted=full[L,c(1,2)]
add=cbind(omitted[,1:2], matrix(data = NA, nrow = nrow(omitted), ncol = ncol(copyframe)-2,))
names(add)[3:ncol(add)]=names(copyframe)[3:ncol(copyframe)]

copyframe=rbind(add,copyframe)
copyframe=copyframe[order(copyframe$CHR,copyframe$POS),]


####adding missing bins to normalized scaled data###
RTs$ID=as.character(paste(RTs$CHR,"_",RTs$POS,sep=""))

L=!full$ID3 %in% RTs$ID
RTs$ID=NULL
omitted=full[L,c(1,2)]
add=cbind(omitted[,1:2], matrix(data = NA, nrow = nrow(omitted), ncol = ncol(RTs)-2,))
names(add)[3:ncol(add)]=names(RTs)[3:ncol(RTs)]

RTs=rbind(add,RTs)
RTs=RTs[order(RTs$CHR,RTs$POS),]


########Binary Heatmap ordered according to rank in S-phase####### (Type "q" to quit)
for(chr in 1:20){

  RTb=subset(copyframe,CHR %in% c(chr))
  RTb.mat=RTb[,as.vector(rank$sample)]
  names(RTb.mat)=as.vector(rank$sample)
  p.binary=heatmap_vd2(RTb.mat,c(1:ncol(RTb.mat)),cols=c("snow2","red"),lim1=0.5, lim2=0.6,col_bias=1,binary=T,legend=T)
  print(p.binary)
  l=readline();if(l=="q") break
}

# ########Normalized scaledd data and binary heatmap with Population data#######

chr=16

RTb=subset(RTs,CHR %in% c(chr))
RTb.mat=RTb[,as.vector(rank$sample)]
names(RTb.mat)=seq(1,nrow(rank),1)
p.trans=heatmap_vd2(RTb.mat,c(1:ncol(RTb.mat)),cols=c("white","snow2","red"),col_bias=2,binary=F,legend=F,lim_buffer=0)


RTb=subset(copyframe,CHR %in% c(chr))
RTb.mat=RTb[,as.vector(rank$sample)]
names(RTb.mat)=seq(1,nrow(rank),1)
p.binary=heatmap_vd2(RTb.mat,c(1:ncol(RTb.mat)),cols=c("snow2","red"),lim1=0.5, lim2=0.6,col_bias=1,binary=T,legend=F)

RTb=subset(RTs,CHR %in% c(chr))
RTb=RTb[,c("CHR","POS","cast.129")]
max.coord=round(max(RTb$POS)/ 20e6,0)
RTplot.bulk=ggplot(RTb, aes(x=POS, y=cast.129,ymin=0,ymax=cast.129)) + geom_ribbon(fill="black") + theme(plot.margin=unit(c(-0.1,0.6,0,-0.3), "cm"))+scale_x_continuous(expand = c(0,0), breaks=seq(0,(max.coord*20e6),20e6),labels=as.character(seq(0,max.coord*20,20)))+geom_hline(yintercept=c(0.05, -0.05), linetype="dashed")

grid.arrange(p.trans,p.binary,RTplot.bulk,layout_matrix=as.matrix(c(1,1,1,2,2,2,3)))

system("mkdir Output_mainscript")

write.table(rank,paste("./Output_mainscript/haploid_rank.txt",sep=""),row.names=F,sep="\t")
write.table(copyframe,paste("./Output_mainscript/haploid_binary.txt",sep=""),row.names=F,sep="\t")
write.table(RTs,paste("./Output_mainscript/haploid_transformed.txt",sep=""),row.names=F,sep="\t")















########################################################################################################
#################################################Remove outlier cells#############################
########################################################################################################

dir=getwd()
dir=paste0(dir,"/Output_mainscript")
setwd(dir)

###########################################################################
##########################Read data#####################################
###########################################################################

copyframe=read.delim("haploid_binary.txt", header=T)
RTbulk=copyframe[,c("CHR","POS","cast.129")]
copyframe$cast.129=NULL

RTs=read.delim("haploid_transformed.txt", header=T)
RTs$cast.129=NULL

#####calculate rank####
rank=data.frame(samples=names(copyframe[,3:ncol(copyframe)]))
rankfun=function(x) length(x[x==1 & is.finite(x)])/length(x[is.finite(x)])*100
rank$rank=apply(copyframe[,3:ncol(copyframe)],2,rankfun)
rank=rank[order(rank$rank),]


########Find outliers#######
hammingdist=function(vec1,vec2) sum(abs(vec1-vec2),na.rm=T)


RTb.mat=copyframe[,as.vector(rank$sample)]

dist.mat=as.data.frame(matrix(NA,ncol(RTb.mat),ncol(RTb.mat)));row.names(dist.mat)=names(RTb.mat);names(dist.mat)=names(RTb.mat)
for(i in 1:ncol(RTb.mat)) {
	for(j in 1:ncol(RTb.mat))
	{
       dist.mat[i,j]=hammingdist(RTb.mat[,i],RTb.mat[,j])
	}
}

dist.mat.p=heatmap_vd2(dist.mat,c(1:ncol(dist.mat)),cols=c("red","orange","yellow"),col_bias=1,binary=F,legend=T)
print(dist.mat.p)

jpeg(file=paste("Hap_outlier_cells_remove.jpeg",sep=""),width=5000,height=5000,res=350)
print(dist.mat.p)
dev.off()



####Remove outliers###

remove=c("G12","D4","E11","G4","C10")

copyframe=copyframe[,-which(names(copyframe) %in% remove)]
RTs=RTs[,-which(names(RTs) %in% remove)]


#####re-calculate rank####
rank=data.frame(samples=names(copyframe[,3:ncol(copyframe)]))
rankfun=function(x) length(x[x==1 & is.finite(x)])/length(x[is.finite(x)])*100
rank$rank=apply(copyframe[,3:ncol(copyframe)],2,rankfun)
rank=rank[order(rank$rank),]




#############FACS accuray calculation#########
rank2=rank[rank$samples %in% c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12"),]
d=density(rank2$rank,adjust=1)
hist(rank2$rank,breaks=seq(0,100,5),freq=T,xlab="Percentage replicated",ylab="Cells",yaxt="n",xaxt="n",border="red",
     xlim=c(0,100),ylim=c(0,12), density=20,col="red")
axis(2, at=seq(0,30,2))
axis(1, at=seq(0,100,10))
par(new=T)
plot(d$y*1.4~d$x,type="l",col="red",xlim=c(0,100),ylim=c(0,0.17),yaxt="n",xaxt="n",ylab="",xlab="",lwd=2)



rank2=rank[!rank$samples %in% c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12","G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12"),]
d=density(rank2$rank,adjust=1)
par(new=T)
hist(rank2$rank,breaks=seq(0,100,5),freq=T,xlab="Percentage replicated",ylab="Cells",yaxt="n",xaxt="n",border="black",
     xlim=c(0,100),ylim=c(0,12), density=10,col="black")
axis(2, at=seq(0,30,2))
axis(1, at=seq(0,100,10))
par(new=T)
plot(d$y/5~d$x,type="l",col="black",xlim=c(0,100),ylim=c(0,0.0066),yaxt="n",xaxt="n",ylab="",xlab="",lwd=2)



rank2=rank[rank$samples %in% c("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12"),]
d=density(rank2$rank,adjust=1)
par(new=T)
hist(rank2$rank,breaks=seq(0,100,5),freq=T,xlab="Percentage replicated",ylab="Cells",yaxt="n",xaxt="n",border="blue",
     xlim=c(0,100),ylim=c(0,12), density=20,col="blue")
axis(2, at=seq(0,30,2))
axis(1, at=seq(0,100,10))
par(new=T)
plot(d$y*1.8~d$x,type="l",col="blue",xlim=c(0,100),ylim=c(0,0.5),yaxt="n",xaxt="n",ylab="",xlab="",lwd=2)




#####write corrected haploid##

copyframe$cast.129=RTbulk$cast.129
RTs$cast.129=RTbulk$cast.129


system("mkdir corrected_files")
write.table(copyframe,paste("./corrected_files/haploid_binary_corrected.txt",sep=""),row.names=F,sep="\t")
write.table(RTs,paste("./corrected_files/haploid_transformed_corrected.txt",sep=""),row.names=F,sep="\t")
write.table(rank,paste("./corrected_files/haploid_rank_corrected.txt",sep=""),row.names=F,sep="\t")




###########MISC###
####heatmap####
for(chr in 1:20){
  cat("\n","chr",chr,"\n")
  RTb=subset(copyframe,CHR %in% c(chr))
  RTb.mat=RTb[,as.vector(rank$sample)]
  names(RTb.mat)=as.vector(rank$sample)
  p.binary=heatmap_vd2(RTb.mat,c(1:ncol(RTb.mat)),cols=c("snow2","red"),lim1=0.5, lim2=0.6,col_bias=1,binary=T,legend=T)
  print(p.binary)
  x=readline() ;if(x=="q") break
}

jpeg(file=paste("Hap_chr1_binary_after.jpeg",sep=""),width=5000,height=3500,res=350)
print(p.binary)
dev.off()


