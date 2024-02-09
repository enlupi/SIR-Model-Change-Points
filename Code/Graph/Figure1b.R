T_max<-10

Delta<-c(1,0,1,0,0,1,0,0,0,0)

set.seed(365)
k<-sum(Delta)
b<-rgamma(k,1,1)
r<-rgamma(k,1,1)

bb<-1/b
rr<-r/(1+r)

png("Change_point.png",600,600,pointsize = 2)
par(mfrow=c(2,1),mar=c(2,2,1,1))
plot(1:T_max,bb[cumsum(Delta)],type="b",pch=17,xlim=c(0,T_max),cex=1.5,lty=1,lwd=2,col="red",xlab = "",ylab = "",xaxt="n",yaxt="n")
mtext("Transmission rate",2,line=0,cex = 2)
# mtext(as.character(1:T_max),1,at=1:T_max)
for(a in which(Delta==1))
{
  abline(v=a,col=1,lty=3,lwd=2)
}
mtext(as.character(1:T_max),1,at=1:T_max,cex = 2,line = 0.5)
mtext(expression(t),1,at=0,cex = 2,line = 0.5)
mtext("=",1,at=0.5,cex = 2,line = 0.5)
par(mar=c(1,2,3,1))
plot(1:T_max,rr[cumsum(Delta)],type="b",pch=19,xlim=c(0,T_max),cex=1.5,lwd=2,lty=1,col="blue",xlab = "",ylab = "",xaxt="n",yaxt="n")
mtext("Removal rate",2,line=0,cex = 2)
for(a in which(Delta==1))
{
  abline(v=a,col=1,lty=3,lwd=2)
}
mtext(as.character(Delta),3,at=1:T_max,line=1.75,cex = 2)
mtext(expression(delta[t]),3,at=0,line=1.75,cex = 2)
mtext("=",3,at=0.5,line=1.75,cex = 2)
mtext(as.character(cumsum(Delta)),3,at=1:T_max,line=0,cex = 2)
mtext(expression(eta[t]),3,at=0,line=0,cex = 2)
mtext("=",3,at=0.5,line=0,cex = 2)
dev.off()