Xuan<-letters[c(1:10)]
Ya<-length(Xuan)

Delta<-rep(c(1,rep(0,24)),4)
Stage<-cumsum(Delta)
K<-sum(Delta)
b<-1/c(0.4,0.4,0.25,0.25)
r<-1/(1-c(0.1,0.25,0.25,0.4))-1
beta<-1/b[Stage]
gamma<-r[Stage]/(1+r[Stage])

Result<-NULL
for(a in Xuan)
{
  Result<-rbind(Result,read.csv(paste(c("Day100III_",a,".csv"),collapse = ""),header = T)[,-1])
}

Stage_dist<-t(apply(Result[1+0:(Ya*10-1)*6,],1,cumsum))
K_max<-max(Stage_dist)
Stage_summary<-NULL
for(k in 1:K_max)
{
  Stage_summary<-rbind(Stage_summary,colMeans(Stage_dist==k))
}
apply(Stage_summary,2,which.max)

b1<-apply(Result[3+0:(Ya*10-1)*6,],2,median)
b2<-apply(Result[4+0:(Ya*10-1)*6,],2,median)
r1<-apply(Result[5+0:(Ya*10-1)*6,],2,median)
r2<-apply(Result[6+0:(Ya*10-1)*6,],2,median)


png("Day100III_Summary.png",400,500)
par(mar=c(4,2,2,1),mfrow=c(2,1))

plot(1:100,(b1),xlim=c(0,100),main = "Transmission rate",ylab="",xlab = "",type="b",pch=18,lty=1,col=2,xaxt="n", yaxt="n")
axis(1,at=0:4*25,labels = as.character(0:4*25),cex.axis=2)
mtext("Day",side=1,line = 3,cex=2)
mtext("Posterior median",side=2,line = 0,cex=2)
for(i in 1+0:3*25)
{
  abline(v=i,col=1,lty=2)
  # lines(i+0:24,beta[i+0:24],col="blue")
}


plot(1:100,(r2),xlim=c(0,100),main = "Removal rate",ylab="",xlab = "",type="b",pch=18,lty=1,col="blue",xaxt="n", yaxt="n")
mtext("Day",side=1,line = 3,cex=2)
mtext("Posterior median",side=2,line = 0,cex=2)
axis(1,at=0:4*25,labels = as.character(0:4*25),cex.axis=2)
for(i in 1+0:3*25)
{
  abline(v=i,col=1,lty=2)
  # lines(i+0:24,gamma[i+0:24],col="blue")
}
dev.off()

P_value_upper_all<-NULL
P_value_lower_all<-NULL
for(i in 0:(Ya*10-1))
{
  Delta_hat<-Result[1+i*6,]
  Delta_pp<-Result[2+i*6,]
  P_value_i_upper<-0*Delta_hat
  P_value_i_lower<-0*Delta_hat
  tau_hat<-which(Delta_hat==1)
  for(j in tau_hat)
  {
    CI_p<-0
    min_a<-1-j
    max_b<-length(Delta_hat)-j
    a<-0
    b<-0
    P_value_i_upper[j+a:b]<-(abs((1-CI_p)+P_value_i_upper[j+a:b])+abs((1-CI_p)-P_value_i_upper[j+a:b]))/2
    CI_p<-sum(Delta_pp[j+a:b])
    P_value_i_lower[j+a:b]<-(abs((1-CI_p)+P_value_i_lower[j+a:b])+abs((1-CI_p)-P_value_i_lower[j+a:b]))/2
    while(CI_p<1)
    {
      pp_zuo<-sum(Delta_pp[j+(a-(a>min_a)):b])
      pp_you<-sum(Delta_pp[j+(a):(b+(b<max_b))])
      if(max(c(pp_zuo,pp_you))>CI_p)
      {
        a<-a-(pp_zuo>pp_you)
        b<-b+(pp_zuo<=pp_you)
        P_value_i_upper[j+a:b]<-(abs((1-CI_p)+P_value_i_upper[j+a:b])+abs((1-CI_p)-P_value_i_upper[j+a:b]))/2
        CI_p<-sum(Delta_pp[j+a:b])
        P_value_i_lower[j+a:b]<-(abs((1-CI_p)+P_value_i_lower[j+a:b])+abs((1-CI_p)-P_value_i_lower[j+a:b]))/2
      }else
      {
        CI_p<-1
      }
    }
  }
  P_value_upper_all<-rbind(P_value_upper_all,(P_value_i_upper))
  P_value_lower_all<-rbind(P_value_lower_all,(P_value_i_lower))
  print(i)
}

png("Day100III_delta_Summary.png",600,400)
par(mar=c(5,5,2,2),mfrow=c(1,1))
plot(2:100,(colMeans(Result[1+0:(Ya*10-1)*6,-1])),xlim=c(0,100),ylim = c(0,1),main = "",ylab="",xlab = "",type="b",pch=18,lwd=2,col=2, xaxs="i", yaxs="i",xaxt="n", yaxt="n",cex=2)
mtext("Day",side=1,line = 3,cex=2)
mtext("Posterior Probability",side=2,line = 3,cex=2)
axis(1,at=0:4*25+c(1,1,1,1,0),labels = as.character(0:4*25+c(1,1,1,1,0)),cex.axis=2)
axis(2,at=0:4/4,labels = as.character(0:4/4),cex.axis=2)
for(i in 1+0:3*25)
{
  abline(v=i,col=1,lty=2)
}
dev.off()

png("Day100III_CP_Summary.png",600,400)
par(mar=c(5,5,2,2),mfrow=c(1,1))
plot(0,0,xlim=c(0,100),ylim = c(0,1),main = "",ylab="",xlab = "",type="n", xaxs="i",xaxt="n", yaxt="n")
mtext("Day",side=1,line = 3,cex=2)
mtext("Coverage Probability",side=2,line = 3,cex=2)
axis(1,at=0:4*25+c(1,1,1,1,0),labels = as.character(0:4*25+c(1,1,1,1,0)),cex.axis=2)
axis(2,at=0:5/5,labels = as.character(0:5/5),cex.axis=2)
for(i in 1+0:3*25)
{
  abline(v=i,col=1,lty=2)
}
alpha_list<-c(0.95,0.9,0.8)
for(alpha in 1:3)
{
  All_CP<-NULL
  for(i in 1:(Ya*10))
  {
    CP<-0*P_value_upper_all[i,]
    CP[P_value_lower_all[i,]>(1-alpha_list[alpha])]<-1
    cross_index<-which((P_value_lower_all[i,]<(1-alpha_list[alpha]))*(P_value_upper_all[i,]>=(1-alpha_list[alpha]))==1)
    CP[cross_index]<-(P_value_upper_all[i,cross_index]-(1-alpha_list[alpha]))/(P_value_upper_all[i,cross_index]-P_value_lower_all[i,cross_index])
    All_CP<-rbind(All_CP,CP)
  }
  lines(2:100,colMeans(All_CP)[-1],col=alpha+1,lty=1,pch=16+alpha,lwd=2)
  abline(h=alpha_list[alpha],col=alpha+1,lty=3,lwd=2)
}

dev.off()
