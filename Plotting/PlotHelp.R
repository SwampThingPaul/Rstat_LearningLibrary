
#Some Example data
library(dataRetrieval)

sdate="1979-10-01"
edate="2010-12-31"
site="02298830"
param="00060"
dat=readNWISdv(site,param,sdate,edate)

#Data manipulation
library(AnalystHelper);# See https://github.com/SwampThingPaul/AnalystHelper

dat$Date=date.fun(dat$Date); #convert "date" field to as.POSIXct
dat$WY=WY(dat$Date,WY.type="Fed");# Determines water year
dat$DOWY=hydro.day(dat$Date,"Fed");# Determines day of water year
dat$discharge_cfs=dat$X_00060_00003
range(dat$WY)

periods=data.frame(WY=seq(1980,2011,1),period=c(rep("1980s",10),rep("1990",10),rep("2000s",10),2010,2011))
dat=merge(dat,periods,"WY")

library(plyr)
period.sum=ddply(subset(dat,period%in%c("1980s","2000s")),c("DOWY","period"),summarise,median.val=median(discharge_cfs,na.rm=T),q10=quantile(discharge_cfs,probs=0.1,na.rm=T),q90=quantile(discharge_cfs,probs=0.9,na.rm=T))

subset(dat,WY==2011)

## Base plotting
ylim.val=c(0,3000);by.y=500;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,366);by.x=90;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)

par(family="serif",mar=c(1.5,3.2,0.5,0.5),oma=c(2,1.5,0.5,0.25),mgp=c(3,1,0));
plot(discharge_cfs~DOWY,dat,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")

with(subset(period.sum,period=="1980s"),shaded.range(DOWY,q10,q90,"red",lty=1))
with(subset(period.sum,period=="1980s"),lines(DOWY,median.val,col="red"))
with(subset(period.sum,period=="2000s"),shaded.range(DOWY,q10,q90,"blue",lty=1))
with(subset(period.sum,period=="2000s"),lines(DOWY,median.val,col="blue"))
with(subset(dat,WY==2010),lines(DOWY,discharge_cfs,lty=1,lwd=3,col="black"))
axis_fun(1,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=2, "Day of Water Year")
mtext(side=2,line=2.75, "Discharge (ft\u00B3 s\u207B\u00B9)")
legend("topleft",legend=c("1980s","2000s","2010"),lty=c(1,1,1),col=c("red","blue","black"),lwd=c(1,1,2),ncol=1,cex=1.25,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
mtext(side=3, "USGS 02298830 (Myakka River, Florida)")
