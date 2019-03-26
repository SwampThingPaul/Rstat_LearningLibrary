## CCA example
## Code was compiled by Paul Julian
## contact infor: pjulian@ufl.edu

#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(vegan)
library(SWMPr)

#Custom functions
axis_fun=function(side,at,at2,labels,cex.axis,line=-0.25,lwd=1){
  axis(side,line=line,at=at,labels=labels,las=1,tcl=-0.6,lty=0,cex.axis=cex.axis);
  axis(side,at=at,labels=F,las=1,tcl=-0.6,lwd=lwd);
  axis(side,at=at2,labels=F,tcl=-0.3,lwd=lwd)
}

##Data
#from vegan package
data(dune)
data(dune.env)


#basic cca see ?cca
## environmental data matrix
mod <- cca(dune ~ A1 + Moisture + Management, dune.env)
plot(mod)

## most info from https://www.fromthebottomoftheheap.net/2012/04/11/customising-vegans-ordination-plots/
with(dune.env,levels(Use))
Use=1:3

scl <- 3 ## scaling = 3
colvec <- c("red2", "green4", "mediumblue")
plot(mod,type="n",scaling=scl)
with(dune.env,points(mod,display="sites",col = colvec[Use],scaling = scl, pch = 21, bg = colvec[Use]))
text(mod, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(dune.env, legend("topright", legend = levels(Use), bty = "n",
                      col = colvec, pch = 21, pt.bg = colvec))


##Customized CCA plot 
##base package ploting only
#scrs=scores(mod)
scrs=scores(mod, display = c("sites", "species"), scaling = scl)
dat=data.frame(scrs$sites)
scrs.bp=scores(mod,display="bp")
mul=ordiArrowMul(scrs.bp)


xlim.val=c(-2.5,2.5);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)#with(scrs, range(species[,1], sites[,1]))
ylim.val=c(-2.5,2.5);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)#with(scrs, range(species[,2], sites[,2]))
plot(CCA2~CCA1,dat,ylim=ylim.val,xlim=xlim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
abline(h=0,v=0,lty=3)
with(dat,points(CCA1,CCA2,pch=21,bg=colvec[Use]));#adds points for scores between CCA 1 and CCA 2
with(scrs,text(species,labels=rownames(species),cex=0.5));#adds points for scores between CCA 1 and CCA 2
arrows(0,0,scrs.bp[,1]*mul,scrs.bp[,2]*mul,length=0.1,angle=15,code=2,col="red");#adds all arrows
text(scrs.bp[,1]*mul,scrs.bp[,2]*mul,labels=rownames(scrs.bp),cex=0.9,col="red")
axis_fun(1,xmaj,xmin,xmaj,1)
axis_fun(2,ymaj,ymin,ymaj,1)
box(lwd=1)
mtext(side=1,line=2.5,"CCA1")
mtext(side=2,line=2.5,"CCA2")






#PCA example
library(SWMPr)
dat=apacpwq
dat=dat[1:500,c("temp","spcond","do_pct","depth","ph","turb")]
dat=na.omit(dat)
pca.all=prcomp(dat,scale=T)
summary(pca.all)
plot(pca.all)
biplot(pca.all)


##PCA plots long hand
eig <- (pca.all$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.pca.all <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)
summary(eig.pca.all)
eig.pca.all

x=barplot(eig.pca.all$variance,ylim=c(0,100),names.arg=1:nrow(eig.pca.all));#each component percent variance
lines(x,eig.pca.all$cumvariance,col="red",lwd=2);#Cumulative variance of all component
barplot(eig.pca.all$eig,names.arg=1:nrow(eig.pca.all));#Eigenvalue for each component

# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
# Variable correlation/coordinates
loadings <- pca.all$rotation
sdev <- pca.all$sdev
var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
head(eig.pca.all)

lab.rownames= rownames(var.coord)

a <- seq(0, 2*pi, length = 100)
plot( cos(a), sin(a), type = 'l', col="gray",ylab=NA,xlab=NA)
abline(h = 0, v = 0, lty = 2)
arrows(0, 0, var.coord[, 1], var.coord[, 2],length = 0.1, angle = 15, code = 2,col="black")
text(1.1*var.coord[,1],1.1*var.coord[,2], labels=lab.rownames, cex = 0.8, adj=1,col="black")
mtext(side=1,line=2,paste0("PCA 1 (",round(eig.pca.all$variance[1],1),"%)"))
mtext(side=2,line=2,paste0("PCA 2 (",round(eig.pca.all$variance[2],1),"%)"))

a <- seq(0, 2*pi, length = 100)
plot( cos(a), sin(a), type = 'l', col="gray",ylab=NA,xlab=NA)
abline(h = 0, v = 0, lty = 2)
arrows(0, 0, var.coord[, 1], var.coord[, 3],length = 0.1, angle = 15, code = 2,col="black")
text(1.1*var.coord[,1],1.1*var.coord[,3], labels=lab.rownames, cex = 0.8, adj=1,col="black")
mtext(side=1,line=2,paste0("PCA 1 (",round(eig.pca.all$variance[1],1),"%)"))
mtext(side=2,line=2,paste0("PCA 3 (",round(eig.pca.all$variance[3],1),"%)"))


#Cos2 plot (quality of representation for variables on the factor map)
var.cos2 <- var.coord^2
#contributions of variables to the PCA
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- data.frame(t(apply(var.cos2,1, contrib, comp.cos2)))
head(var.contrib[, 1:4])

var.contrib.PCA1=var.contrib[order(-var.contrib$PC1),]

var.contrib.PCA1.2=var.contrib[,c("PC1","PC2")]
var.contrib.PCA1.2$PerComp=with(var.contrib.PCA1.2,(PC1*eig.pca.all$eig[1])+(PC2*eig.pca.all$eig[2]))
var.contrib.PCA1.2=var.contrib.PCA1.2[order(-var.contrib.PCA1.2$PerComp),]

var.contrib.PCA1.3=var.contrib[,c("PC1","PC3")]
var.contrib.PCA1.3$PerComp=with(var.contrib.PCA1.3,(PC1*eig.pca.all$eig[1])+(PC3*eig.pca.all$eig[3]))
var.contrib.PCA1.3=var.contrib.PCA1.3[order(-var.contrib.PCA1.3$PerComp),]

barplot(var.contrib.PCA1.2$PerComp,names.arg=rownames(var.contrib.PCA1.2))
cut.off=(length(rownames(var.contrib))*eig.pca.all$eig[1])+(length(rownames(var.contrib))*eig.pca.all$eig[2])
abline(h=cut.off,col="red",lty=2,lwd=2,xpd=F)

barplot(var.contrib.PCA1.3$PerComp,names.arg=rownames(var.contrib.PCA1.3))
cut.off=(length(rownames(var.contrib))*eig.pca.all$eig[1])+(length(rownames(var.contrib))*eig.pca.all$eig[2])
abline(h=cut.off,col="red",lty=2,lwd=2,xpd=F)
