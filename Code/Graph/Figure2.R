T_max<-10

Delta<-c(1,0,1,0,0,1,0,0,0,0)
Delta<-cbind(Delta,c(1,0,1,0,0,1,0,0,1,0))
Delta<-cbind(Delta,c(1,0,0,0,0,1,0,0,0,0))
Delta<-cbind(Delta,c(1,0,1,0,0,0,1,0,0,0))
kk<-c("Initialization","Add","Delete","Swap")

library(shape)
png("Add_drop_swap.png",600,600)
par(mfrow=c(4,1),mar=c(1,2,1,1))
for(al in 1:4)
{
  plot(0,0,type="n",xlim=c(0,T_max+1),ylim = c(-1.2,1),lty=1,xlab = "",ylab = "",xaxt="n",yaxt="n",axes = F)
  Arrows(0,0,T_max+1,0,lwd=2)
  # text(0,T_max+1,'t',cex=3)
  text((1+T_max)/2,0.7,kk[al],adj=c(0.5,0),col=1,cex=3)
  text(0,-0.1,expression(t),adj=c(0.5,1),col=1,cex=3)
  text(0,-0.6,expression(delta[t]^"*"),adj=c(0.5,1),col=1,cex=3)
  for(i in 1:T_max)
  {
    
    if(Delta[i,al]==0)
    {
      segments(i,0,i,0.3,lwd=2,col=1)
      text(i,-0.1,as.character(i),adj=c(0.5,1),col=1,cex=3)
      text(i,-0.7,as.character(Delta[i,al]),adj=c(0.5,1),col=1,cex=3)
    }
    else
    {
      segments(i,0,i,0.3,lwd=2,col=2)
      points(i,0.3,cex=2,pch=17,col=2)
      text(i,-0.1,as.character(i),adj=c(0.5,1),col=2,cex=3)
      text(i,-0.7,as.character(Delta[i,al]),adj=c(0.5,1),col=2,cex=3)
    }
  }
}
dev.off()