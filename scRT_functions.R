
####Author: Vishnu Dileep
####Updated: Jan 2018
###########################################################################
##########################Single-cell RT Functions#####################################
###########################################################################
rescaleDatasets = function(RT) {
	for (i in 3:ncol(RT)) {						# 4:ncol(RT) if probeID present
		RT[,i] = RT[,i] * 1.59 / IQR(RT[,i],na.rm=T)	# Scale to BG01ESC.R1 IQR
		cat("Scaling dataset", i-2, "/", ncol(RT)-2, "\n")
	}
	return(RT)
}    




heatmap_vd2= function (mat,plot_incl,cols=c("green","black","red"),col_bias=1,legend_range=0.01,lim1=NULL,lim2=NULL,binary=F,legend=F,lim_buffer=0){
	
	###transform matrix for levelplot####
	mat=as.matrix(mat[,plot_incl])
    mat <- mat[,seq(from=ncol(mat),to=1,by=-1)]
    
    ###setcolors####
	if(is.null(lim1)){
	lim1=quantile(stack(as.data.frame(mat[,plot_incl]))$values,legend_range/100,na.rm=T)-lim_buffer
    lim2=quantile(stack(as.data.frame(mat[,plot_incl]))$values,1-(legend_range/100),na.rm=T)+lim_buffer}
    	
    mycols=colorRampPalette(cols,col_bias)( length(seq(lim1,lim2,(lim2-lim1)/1000) ))
    
    ###set axis##
    max.coord=round(nrow(mat)*50000/ 20e6,0)
    
    ###set legend param##
    if(legend==F) colkey_param=F
    if(legend==T) colkey_param=list(at=seq(lim1,lim2,(lim2-lim1)/1000), labels=list(at=round(seq(lim1,lim2,(lim2-lim1)/10),1)), height=0.5)
    
      
  ###Plot###  
if(binary==T){
    p=levelplot(mat, xlab="", ylab="", col.regions=mycols, cuts=2, at=seq(0,1,0.5),aspect="fill",useRaster=T,
                colorkey=colkey_param, 
                scales=list(x=list(at=seq(0,(max.coord*20e6)/50000,20e6/50000),labels=as.character(seq(0,max.coord*20,20)))),
                par.settings = list(panel.background=list(col="white"),layout.heights = list( 
                        top.padding = -2, 
                        main.key.padding = 0, 
                        key.axis.padding = 0, 
                        axis.xlab.padding = 0, 
                        xlab.key.padding = 0, 
                        key.sub.padding = 0,
                        bottom.padding=0 ), 
                layout.widths = list( 
                        left.padding = 0, 
                        key.ylab.padding = 0, 
                        ylab.axis.padding = 0, 
                        axis.key.padding = 0, 
                        right.padding = 0) 
                )  )
           }
        

 if(binary==F){
    p=levelplot(mat, xlab="", ylab="", col.regions=mycols, at=seq(lim1,lim2,(lim2-lim1)/1000),aspect="fill",useRaster=T,
                colorkey=colkey_param, 
                scales=list(x=list(at=seq(0,(max.coord*20e6)/50000,20e6/50000), labels=as.character(seq(0,max.coord*20,20)))), 
                par.settings=list(panel.background=list(col="white"),layout.heights = list( 
                        top.padding = -2, 
                        main.key.padding = 0, 
                        key.axis.padding = 0, 
                        axis.xlab.padding = 0, 
                        xlab.key.padding = 0, 
                        key.sub.padding = 0,
                        bottom.padding=0 ), 
                layout.widths = list( 
                        left.padding = 0, 
                        key.ylab.padding = 0, 
                        ylab.axis.padding = 0, 
                        axis.key.padding = 0, 
                        right.padding = 0)))         
               }
    return(p)        
	}
	
	
	
	####check segmentation###

check.seg= function(RT,segframe,binary=F,samples=NULL,chrs=1,x1=2e6,x2=99e6){
if(is.null(samples)) samples=names(RT)[-1*1:2]
             ###plot##
 m=as.matrix(rbind(c(1,1,1),c(1,1,1),c(2,2,2)))
 layout(m)
 #layout.show(3)
par(mar=c(2,3,1,0.5),ask=T, omi = c(0.4,0,0.1,0.5))


for(sample in samples){
for(chr in chrs){
seg= segframe[segframe$CHR==chr, c(1,2,which(names(segframe)==sample))]

data=subset(RT,CHR==chr)
data2=subset(RTbulk,CHR==chr)

data=merge(data,data2,by=c("CHR","POS"))
data=data[order(data$CHR,data$POS),]


                      ## Plot single cell data###
y1=quantile(data[,which(names(data)==sample)],0.001,na.rm=T)
y2=quantile(data[,which(names(data)==sample)],0.999,na.rm=T)
plot(data[,which(names(data)==sample)]~data$POS,axes=FALSE,xlim=c(x1,x2), ylim=c(y1,y2), xlab="", ylab="",pch=19,type="p",col="black",cex=0.3, main=paste("SC RT Copy number ",sample,"CHR",chr,sep=" "))
axis(2,col="black",las=1)  ## las=1 makes horizontal labels
mtext("rpm",side=2,line=2)
box()

## Allow a second plot on the same graph
par(new=TRUE)


## Plot the second plot and put axis scale on right
if(binary==F) plot(seg[,3]~seg$POS, col="orangered",ylim=c(y1,y2), xlab="", ylab="",axes=FALSE, type="l",xlim=c(x1,x2),lwd=1,lty=1)
if(binary==T) plot(seg[,3]~seg$POS, col="orangered",ylim=c(-0.2,1.2), xlab="", ylab="",axes=FALSE, type="l",xlim=c(x1,x2),lwd=1,lty=1)
## a little farther out (line=4) to make room for labels
mtext("Copy number",side=4,col="orangered",line=2) 
axis(4, col="red",col.axis="red",las=1)

## Draw the time axis
axis(1,at=seq(0,round(max(data$POS)/1e6,0),20)*1e6,labels=as.character(seq(0,round(max(data$POS)/1e6,0),20)))

                             ######Bulk plot###
plot(data2$cast.129~data2$POS,axes=FALSE,xlim=c(x1,x2), ylim=c(-4.5,4.5), xlab="", ylab="",pch=19,type="h",col="black",cex=0.3, main=paste("Bulk RT ",sample,"CHR",chr,sep=" "))
axis(2,col="black",las=1)  ## las=1 makes horizontal labels
mtext("RT",side=2,line=2)
box()

## Allow a second plot on the same graph
par(new=TRUE)

if(binary==F) plot(seg[,3]~seg$POS, col="orangered",ylim=c(y1,y2), xlab="", ylab="",axes=FALSE, type="l",xlim=c(x1,x2),lwd=1,lty=1)
if(binary==T) plot(seg[,3]~seg$POS, col="orangered",ylim=c(-0.2,1.2), xlab="", ylab="",axes=FALSE, type="l",xlim=c(x1,x2),lwd=1,lty=1)
## a little farther out (line=4) to make room for labels
mtext("Copy number",side=4,col="orangered",line=2) 
axis(4, col="red",col.axis="red",las=1)

## Draw the time axis
axis(1,at=seq(0,round(max(data$POS)/1e6,0),20)*1e6,labels=as.character(seq(0,round(max(data$POS)/1e6,0),20)))

 }
 }
}




