I_all<-NULL
R_all<-NULL
S_all<-NULL
NI_all<-NULL
NR_all<-NULL

for(iter in 1:100)
{
  N<-1000000
  T_max<-100
  
  Delta<-rep(c(1,rep(0,24)),4)
  Stage<-cumsum(Delta)
  K<-sum(Delta)
  
  b<-1/c(0.4,0.4,0.25,0.25)
  r<-1/(1-c(0.1,0.25,0.25,0.4))-1
  
  beta<-1/b[Stage]
  gamma<-r[Stage]/(1+r[Stage])
  
  I_current<-50
  R_current<-0
  S_current<-N-I_current-R_current
  
  N_infect_obs<-I_current
  N_recovery_obs<-R_current
  
  for(t in 1:T_max)
  {
    PI_current<-I_current/N
    p_beta<-1-exp(-beta[t]*PI_current)
    New_infect<-rbinom(1,S_current,p_beta)
    New_recovery<-rbinom(1,I_current,gamma[t])
    
    N_infect_obs<-c(N_infect_obs,New_infect)
    N_recovery_obs<-c(N_recovery_obs,New_recovery)
    S_current<-S_current-New_infect
    I_current<-I_current-New_recovery+New_infect
    R_current<-R_current+New_recovery
  }
  
  S_obs<-c(N-cumsum(N_infect_obs))[1+0:T_max]
  I_obs<-c(cumsum(N_infect_obs-N_recovery_obs))[1+0:T_max]
  R_obs<-c(cumsum(N_recovery_obs))[1+0:T_max]
  PI_obs<-I_obs/N
  N_infect_obs<-N_infect_obs[-1]
  N_recovery_obs<-N_recovery_obs[-1]
  
  I_all<-rbind(I_all,I_obs)
  R_all<-rbind(R_all,R_obs)
  S_all<-rbind(S_all,S_obs)
  NI_all<-rbind(NI_all,N_infect_obs)
  NR_all<-rbind(NR_all,N_recovery_obs)
}

png("Mean_II.png",600,600)
par(mar=c(5,5,1,2))
plot(0,0,type = "n",xlim = c(-0.5,T_max+0.5),ylim = c(0,1),xlab="",ylab="",axes = T, xaxs="i", yaxs="i",xaxt="n", yaxt="n")

mtext("Day",side=1,line = 3,cex=2)
mtext("Proportion",side=2,line = 3,cex=2)
rect(-0.5,0,T_max+0.5,1,lwd=3)
axis(1,at=0:4*25,labels = as.character(0:4*25),cex.axis=2)
axis(2,at=0:5/5,labels = as.character(0:5/5),cex.axis=2)
for(iter in 1:100)
{
  for(t in 0:T_max)
  {
    rect(t-0.5,0,t+0.5,S_all[iter,t+1]/N,col=rgb(0,1,0,alpha = 0.01,maxColorValue = 1),border = rgb(0,1,0,alpha = 0.01,maxColorValue = 1))
    rect(t-0.5,S_all[iter,t+1]/N,t+0.5,(S_all[iter,t+1]+I_all[iter,t+1])/N,col=rgb(1,0,0,alpha = 0.01,maxColorValue = 1),border = rgb(1,0,0,alpha = 0.01,maxColorValue = 1))
    rect(t-0.5,(S_all[iter,t+1]+I_all[iter,t+1])/N,t+0.5,1,col=rgb(0,0,1,alpha = 0.01,maxColorValue = 1),border = rgb(0,0,1,alpha = 0.01,maxColorValue = 1))
  }
}
legend("topright",col=c("green","red","blue"),legend = c("S","I","R"),pch=15,cex=2)
dev.off()

####################################################


I_all<-NULL
R_all<-NULL
S_all<-NULL
NI_all<-NULL
NR_all<-NULL

for(iter in 1:100)
{
  N<-1000000
  T_max<-100
  
  Delta<-rep(c(1,rep(0,24)),4)
  Stage<-cumsum(Delta)
  K<-sum(Delta)
  
  b<-1/c(0.5,0.3,0.4,0.2)
  r<-1/(1-c(0.1,0.3,0.2,0.4))-1
  
  beta<-1/b[Stage]
  gamma<-r[Stage]/(1+r[Stage])
  
  I_current<-50
  R_current<-0
  S_current<-N-I_current-R_current
  
  N_infect_obs<-I_current
  N_recovery_obs<-R_current
  
  for(t in 1:T_max)
  {
    PI_current<-I_current/N
    p_beta<-1-exp(-beta[t]*PI_current)
    New_infect<-rbinom(1,S_current,p_beta)
    New_recovery<-rbinom(1,I_current,gamma[t])
    
    N_infect_obs<-c(N_infect_obs,New_infect)
    N_recovery_obs<-c(N_recovery_obs,New_recovery)
    S_current<-S_current-New_infect
    I_current<-I_current-New_recovery+New_infect
    R_current<-R_current+New_recovery
  }
  
  S_obs<-c(N-cumsum(N_infect_obs))[1+0:T_max]
  I_obs<-c(cumsum(N_infect_obs-N_recovery_obs))[1+0:T_max]
  R_obs<-c(cumsum(N_recovery_obs))[1+0:T_max]
  PI_obs<-I_obs/N
  N_infect_obs<-N_infect_obs[-1]
  N_recovery_obs<-N_recovery_obs[-1]
  
  I_all<-rbind(I_all,I_obs)
  R_all<-rbind(R_all,R_obs)
  S_all<-rbind(S_all,S_obs)
  NI_all<-rbind(NI_all,N_infect_obs)
  NR_all<-rbind(NR_all,N_recovery_obs)
}

png("Mean_III.png",600,600)
par(mar=c(5,5,1,2))
plot(0,0,type = "n",xlim = c(-0.5,T_max+0.5),ylim = c(0,1),xlab="",ylab="",axes = T, xaxs="i", yaxs="i",xaxt="n", yaxt="n")

mtext("Day",side=1,line = 3,cex=2)
mtext("Proportion",side=2,line = 3,cex=2)
rect(-0.5,0,T_max+0.5,1,lwd=3)
axis(1,at=0:4*25,labels = as.character(0:4*25),cex.axis=2)
axis(2,at=0:5/5,labels = as.character(0:5/5),cex.axis=2)
for(iter in 1:100)
{
  for(t in 0:T_max)
  {
    rect(t-0.5,0,t+0.5,S_all[iter,t+1]/N,col=rgb(0,1,0,alpha = 0.01,maxColorValue = 1),border = rgb(0,1,0,alpha = 0.01,maxColorValue = 1))
    rect(t-0.5,S_all[iter,t+1]/N,t+0.5,(S_all[iter,t+1]+I_all[iter,t+1])/N,col=rgb(1,0,0,alpha = 0.01,maxColorValue = 1),border = rgb(1,0,0,alpha = 0.01,maxColorValue = 1))
    rect(t-0.5,(S_all[iter,t+1]+I_all[iter,t+1])/N,t+0.5,1,col=rgb(0,0,1,alpha = 0.01,maxColorValue = 1),border = rgb(0,0,1,alpha = 0.01,maxColorValue = 1))
  }
}
legend("topright",col=c("green","red","blue"),legend = c("S","I","R"),pch=15,cex=2)
dev.off()


#############################################

I_all<-NULL
R_all<-NULL
S_all<-NULL
NI_all<-NULL
NR_all<-NULL

for(iter in 1:100)
{
  N<-1000000
  T_max<-100
  
  Delta<-rep(c(1,rep(0,24)),4)
  Stage<-cumsum(Delta)
  K<-sum(Delta)
  
  b<-1/c(0.3,0.4,0.25,0.2)
  r<-1/(1-c(0.05,0.15,0.2,0.25))-1
  
  beta<-1/b[Stage]
  gamma<-r[Stage]/(1+r[Stage])
  
  I_current<-50
  R_current<-0
  S_current<-N-I_current-R_current
  
  N_infect_obs<-I_current
  N_recovery_obs<-R_current
  
  for(t in 1:T_max)
  {
    PI_current<-I_current/N
    p_beta<-1-exp(-beta[t]*PI_current)
    New_infect<-rbinom(1,S_current,p_beta)
    New_recovery<-rbinom(1,I_current,gamma[t])
    
    N_infect_obs<-c(N_infect_obs,New_infect)
    N_recovery_obs<-c(N_recovery_obs,New_recovery)
    S_current<-S_current-New_infect
    I_current<-I_current-New_recovery+New_infect
    R_current<-R_current+New_recovery
  }
  
  S_obs<-c(N-cumsum(N_infect_obs))[1+0:T_max]
  I_obs<-c(cumsum(N_infect_obs-N_recovery_obs))[1+0:T_max]
  R_obs<-c(cumsum(N_recovery_obs))[1+0:T_max]
  PI_obs<-I_obs/N
  N_infect_obs<-N_infect_obs[-1]
  N_recovery_obs<-N_recovery_obs[-1]
  
  I_all<-rbind(I_all,I_obs)
  R_all<-rbind(R_all,R_obs)
  S_all<-rbind(S_all,S_obs)
  NI_all<-rbind(NI_all,N_infect_obs)
  NR_all<-rbind(NR_all,N_recovery_obs)
}

png("Mean_I.png",600,600)
par(mar=c(5,5,1,2))
plot(0,0,type = "n",xlim = c(-0.5,T_max+0.5),ylim = c(0,1),xlab="",ylab="",axes = T, xaxs="i", yaxs="i",xaxt="n", yaxt="n")

mtext("Day",side=1,line = 3,cex=2)
mtext("Proportion",side=2,line = 3,cex=2)
rect(-0.5,0,T_max+0.5,1,lwd=3)
axis(1,at=0:4*25,labels = as.character(0:4*25),cex.axis=2)
axis(2,at=0:5/5,labels = as.character(0:5/5),cex.axis=2)
for(iter in 1:100)
{
  for(t in 0:T_max)
  {
    rect(t-0.5,0,t+0.5,S_all[iter,t+1]/N,col=rgb(0,1,0,alpha = 0.01,maxColorValue = 1),border = rgb(0,1,0,alpha = 0.01,maxColorValue = 1))
    rect(t-0.5,S_all[iter,t+1]/N,t+0.5,(S_all[iter,t+1]+I_all[iter,t+1])/N,col=rgb(1,0,0,alpha = 0.01,maxColorValue = 1),border = rgb(1,0,0,alpha = 0.01,maxColorValue = 1))
    rect(t-0.5,(S_all[iter,t+1]+I_all[iter,t+1])/N,t+0.5,1,col=rgb(0,0,1,alpha = 0.01,maxColorValue = 1),border = rgb(0,0,1,alpha = 0.01,maxColorValue = 1))
  }
}
legend("topright",col=c("green","red","blue"),legend = c("S","I","R"),pch=15,cex=2)
dev.off()

