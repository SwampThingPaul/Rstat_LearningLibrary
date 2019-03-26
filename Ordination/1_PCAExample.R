
## Example of PCA/RDA code
##
## Code was compiled by Paul Julian
## contact infor: pjulian@ufl.edu

#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(vegan)

#Custom functions
# Makes pretty axes when using base R-plotting
axis_fun=function(side,at,at2,labels,cex.axis,line=-0.25,lwd=1){
  axis(side,line=line,at=at,labels=labels,las=1,tcl=-0.6,lty=0,cex.axis=cex.axis);
  axis(side,at=at,labels=F,las=1,tcl=-0.6,lwd=lwd);
  axis(side,at=at2,labels=F,tcl=-0.3,lwd=lwd)
}

#Import data (using stock data for example)
data(dune)
data(dune.env)

## PCA Analysis
my.rda=rda(dune,scale=T);#  if scale=T then the analysis is identifical to PCA

biplot(my.rda);#stock biplot

# This pulls all of the eigenvalue and 
# estimates the variance for each component
eig <- my.rda$CA$eig
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.pca <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)

#Screeplots
par(family="serif",mar=c(0.5,1.5,1,1),oma=c(3.5,2,0.25,0.25));

ylim.val=c(0,105);by.y=25;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2);#set y limit and delineates the major and minor ticks
x=barplot(eig.pca$variance,ylim=ylim.val,col="white",border=0,yaxt="n")# inital plot to get the measurements
abline(h=ymaj,lty=3,col="grey")#makes vertical lines from y axis
x=barplot(eig.pca$variance,ylim=ylim.val,col="grey",yaxt="n",add=T)# the real plot that matters
lines(x,eig.pca$cumvariance,col="indianred1",lwd=2)# adds the cumulative variance for each factor
points(x,eig.pca$cumvariance,pch=21,bg="indianred1",cex=1.25)
axis_fun(1,x,x,seq(1,length(x),1),1)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Principal Components")
mtext(side=2,line=2,"Percentage of Variances")
legend.text=c("Absolute","Cumulative");#helper vaiable for legend
pt.col=c("grey","indianred1")#helper vaiable for legend
legend("topleft",legend=legend.text,pch=c(22,21),pt.bg=pt.col,col=c("black",pt.col[2]),lty=c(0,1),lwd=1.5,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend("topleft",legend=legend.text,pch=c(22,21),pt.bg=pt.col,col="black",lty=0,lwd=0.5,pt.cex=1.55,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)

ylim.val=c(0,8);by.y=2;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2)
x=barplot(eig.pca$eig,ylim=ylim.val,col="grey",yaxt="n")
abline(h=ymaj,lty=3,col="grey")
x=barplot(eig.pca$eig,ylim=ylim.val,col="grey",yaxt="n",add=T)
axis_fun(1,line=-0.7,x,x,seq(1,length(x),1),0.7)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Principal Components")
mtext(side=2,line=1.5,"Eigenvalue")


#biplot
scrs=scores(my.rda,display=c("sites","species"));#extracts scores from PCA analysis

#Basic biplot you get from biplot()
xlim.val=c(-3,3);by.x=1;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);# xlim., major and minor ticks all in one line
ylim.val=c(-3,2);by.y=1;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);# ylim., major and minor ticks all in one line
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);#essentially an empty plot
abline(h=0,v=0,lty=3,col="grey");# makes the fancy 0 axis lines
points(scrs$sites,pch=21,bg="grey",cex=1,lwd=0.5); #plots the points
arrows(0,0,scrs$species[,1],scrs$species[,2],length = 0.05, angle = 15, code = 2,col="indianred1",lwd=1.5);# makes the arrows
with(scrs,text(species[,1]+0.25,species[,2],labels=rownames(species),cex=0.5));#adds labels to the arrows; 
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); #adds x axis ticks
axis_fun(2,ymaj,ymin,format(ymaj),1); #adds y axis ticks
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2,paste0("PCA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance




# if you prefer to add some additional information ordiellipse 
# is a great way to highlight different attributes
# Can also code the points to difference values for the same effect. 

xlim.val=c(-3,3);by.x=1;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(-3,2);by.y=1;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
abline(h=0,v=0,lty=3,col="grey")
points(scrs$sites,pch=21,bg="grey",cex=1,lwd=0.5)
x=ordiellipse(my.rda,group=dune.env$Management,draw="polygon",label=F,col=c("dodgerblue1","indianred1","forestgreen","goldenrod"),border=T,cex=0.8)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1)
axis_fun(2,ymaj,ymin,format(ymaj),1)
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2,paste0("PCA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance


# Code the points to be different shape based on management field
dune.env.points=merge(data.frame(id=1:nrow(dune.env),Management=dune.env$Management),
                      data.frame(Management=c("SF","BF","HF","NM"),
                                 point.num=c(21,22,23,24),
                                 point.col=c("dodgerblue1","indianred1","forestgreen","goldenrod")))

# you have to trick the data set because when you merge the management data
# it reorders the data which is why I made the "id" variable. 
dune.env.points=dune.env.points[order(dune.env.points$id),];

# Biplot with ellispes 
xlim.val=c(-3,3);by.x=1;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(-3,2);by.y=1;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
abline(h=0,v=0,lty=3,col="grey")
points(scrs$sites,pch=dune.env.points$point.num,bg="grey",cex=1,lwd=0.5)
x=ordiellipse(my.rda,group=dune.env$Management,draw="polygon",label=F,col=c("dodgerblue1","indianred1","forestgreen","goldenrod"),border=T,cex=0.8)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1)
axis_fun(2,ymaj,ymin,format(ymaj),1)
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2,paste0("PCA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance
legend.text=c("SF","BF","HF","NM");#helper vaiable for legend
legend("topleft",legend=legend.text,pch=c(21,22,23,24),pt.bg="grey",pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)


# Biplot different shape and colors
xlim.val=c(-3,3);by.x=1;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(-3,2);by.y=1;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
abline(h=0,v=0,lty=3,col="grey")
points(scrs$sites,pch=dune.env.points$point.num,bg=adjustcolor(dune.env.points$point.col,0.5),cex=1,lwd=0.5);#adjusts color and shape
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1)
axis_fun(2,ymaj,ymin,format(ymaj),1)
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2,paste0("PCA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance
legend.text=c("SF","BF","HF","NM");#helper vaiable for legend
pt.col=adjustcolor(c("dodgerblue1","indianred1","forestgreen","goldenrod"),0.5)
legend("topleft",legend=legend.text,pch=c(21,22,23,24),pt.bg=pt.col,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
