library(shape)
png("SIR_graph.png",450,300)
par(mfrow=c(1,1),mar=c(1,1,1,1))
plot(0,0,type="n",xlim=c(0,7.5),ylim=c(c(-2.5,2.5)),axes = F)
for(i in 0:2)
{
  rect(0+2.5*i,1.5,1+2.5*i,2.5,lwd=2)
  text(0.5+2.5*i,2,parse(text=paste("S[",as.character(i),"]")),cex=2.5)
  
  rect(0+2.5*i,-0.5,1+2.5*i,0.5,border = "red",lwd=2)
  text(0.5+2.5*i,0,parse(text=paste("I[",as.character(i),"]")),col="red",cex=2.5)
  
  rect(0+2.5*i,-1.5,1+2.5*i,-2.5,border = "blue",lwd=2)
  text(0.5+2.5*i,-2,parse(text=paste("R[",as.character(i),"]")),col="blue",cex=2.5)
}

for(i in 0:1)
{
  Arrows(1+2.5*i,2,2.35+2.5*i,2,lwd=2)
  text(1.75+2.5*i,2,parse(text=paste("S[",as.character(i),"]-Delta")),adj=c(0.8,-0.2),cex=1.5)
  text(1.75+2.5*i,2,parse(text=paste("I[",as.character(i+1),"]")),adj=c(-0.8,-0.2),cex=1.5)
  
  Arrows(1+2.5*i,1.5,2.35+2.5*i,0.6,col="red",lwd=2)
  text(1.75+2.5*i,1+0.05*2,parse(text=paste("Delta")),adj=c(0.8,-0.2),srt=0,col="red",cex=1.5)
  text(1.75+2.5*i,1,parse(text=paste("I[",as.character(i+1),"]")),adj=c(-0.2,-0.2),srt=0,col="red",cex=1.5)
  
  Arrows(1+2.5*i,0,2.35+2.5*i,0,col="red",lwd=2)
  text(1.75+2.5*i,0,parse(text=paste("I[",as.character(i),"]-Delta")),adj=c(0.8,-0.2),col="red",cex=1.5)
  text(1.75+2.5*i,0,parse(text=paste("R[",as.character(i+1),"]")),adj=c(-0.4,-0.2),col="red",cex=1.5)
  
  Arrows(1+2.5*i,-0.5,2.35+2.5*i,-1.4,col="blue",lwd=2)
  text(1.75+2.5*i,-1+0.05*2,parse(text=paste("Delta")),adj=c(0.8,-0.2),srt=0,col="blue",cex=1.5)
  text(1.75+2.5*i,-1,parse(text=paste("R[",as.character(i+1),"]")),adj=c(-0.2,-0.2),srt=0,col="blue",cex=1.5)
  
  Arrows(1+2.5*i,-2,2.35+2.5*i,-2,col="blue",lwd=2,cex=1.5)
  text(1.75+2.5*i,-2,parse(text=paste("R[",as.character(i),"]")),adj=c(0.5,-0.2),col="blue",cex=1.5)
}

for(i in 15:17/18*7.5)
{
  points(i,2,pch=19,cex=1.2)
  points(i,0,pch=19,cex=1.2,col="red")
  points(i,-2,pch=19,cex=1.2,col="blue")
}
dev.off()