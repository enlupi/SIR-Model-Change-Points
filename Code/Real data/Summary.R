X<-read.csv("Singapore_new.csv",header=T)
Active<-cumsum(X$Confirmed-X$Recovered)
wa<-min(which(Active[-length(Active)]>X$Recovered[-1]))
N_infect<-X$Confirmed
N_infect<-smooth(N_infect)
N_recovery<-c(rep(0,wa),X$Recovered[-c(1:wa)])
N_recovery<-smooth(N_recovery)

N<-5930134
T_max<-length(N_infect)
T_pred<-0

Burnin<-5000
AAA<-1000
thining<-10

N_infect_obs<-N_infect[1:T_max]
N_recovery_obs<-N_recovery[1:T_max]
S_obs<-c(N-cumsum(N_infect))[1:T_max]
I_obs<-c(cumsum(N_infect-N_recovery))[1:T_max]
R_obs<-c(cumsum(N_recovery))[1:T_max]
PI_obs<-I_obs/N
N_infect_obs<-N_infect_obs[-1]
N_recovery_obs<-N_recovery_obs[-1]
T_max<-T_max-1

############################################################
Beta_all<-read.csv("Beta_all.csv")
B_all<-read.csv("B_all.csv")
Gamma_all<-read.csv("Gamma_all.csv")
R_all<-read.csv("R_all.csv")
Delta_all<-read.csv("Delta_all.csv")

B_median<-apply(B_all,2,median)
B_mean<-apply(B_all,2,mean)
B_sd<-apply(B_all,2,sd)
B_CI<-apply(B_all,2,quantile,prob=c(0.025,0.975))

R_median<-apply(R_all,2,median)
R_mean<-apply(R_all,2,mean)
R_sd<-apply(R_all,2,sd)
R_CI<-apply(R_all,2,quantile,prob=c(0.025,0.975))



Stage_all<-(t(apply(Delta_all,1,cumsum)))
PPM<-matrix(0,ncol = ncol(Stage_all),nrow = ncol(Stage_all))
for(i in 1:nrow(Stage_all))
{
  PPM<-PPM+(Stage_all[i,row(PPM)]==Stage_all[i,col(PPM)])/nrow(Stage_all)
}
Delta_final<-rep(0,T_max)
Delta_final[1]<-1
Stage_hat<-cumsum(Delta_final)
Continue_index_add_drop<-TRUE
Continue_index_swap<-TRUE
Current_loss<-sum(abs((Stage_hat[row(PPM)]==Stage_hat[col(PPM)])-PPM))
all_candidate_index<-2:T_max
while ((Continue_index_add_drop+Continue_index_swap)>0) 
{
  All_loss<-NULL
  all_candidate_index<-2:T_max
  for(i in all_candidate_index)
  {
    Delta_candidate<-Delta_final
    Delta_candidate[i]<-1-Delta_candidate[i]
    
    Stage_candidate<-cumsum(Delta_candidate)
    Candidate_loss<-sum(abs((Stage_candidate[row(PPM)]==Stage_candidate[col(PPM)])-PPM))
    All_loss<-c(All_loss,Candidate_loss)
  }
  if(min(All_loss)<Current_loss)
  {
    Current_loss<-min(All_loss)
    Delta_final[all_candidate_index[which.min(All_loss)]]<-1-Delta_final[all_candidate_index[which.min(All_loss)]]
    Continue_index_add_drop<-TRUE
  }else
  {
    Continue_index_add_drop<-FALSE
  }
  
  if(sum(Delta_final)%in%(2:(T_max-1)))
  {
    All_loss<-NULL
    all_candidate_index<-which((Delta_final[(2:(T_max-1))]-Delta_final[(2:(T_max-1))+1])!=0)+1
    for(i in all_candidate_index)
    {
      Delta_candidate<-Delta_final
      Delta_candidate[i+0:1]<-1-Delta_candidate[i+0:1]
      
      Stage_candidate<-cumsum(Delta_candidate)
      Candidate_loss<-sum(abs((Stage_candidate[row(PPM)]==Stage_candidate[col(PPM)])-PPM))
      All_loss<-c(All_loss,Candidate_loss)
    }
    if(min(All_loss)<Current_loss)
    {
      Current_loss<-min(All_loss)
      Delta_final[all_candidate_index[which.min(All_loss)]+0:1]<-1-Delta_final[all_candidate_index[which.min(All_loss)]+0:1]
      Continue_index_swap<-TRUE
    }
    else
    {
      Continue_index_swap<-FALSE
    }
  }
}
Delta_mean<-colMeans(Delta_all)
# cbind(Delta_final,Delta_mean)




Delta_hat<-Delta_final
Delta_pp<-Delta_mean
P_value_i_upper<-0*Delta_hat
P_value_i_lower<-0*Delta_hat
tau_hat<-which(Delta_hat==1)
CI_tau<-list()
for(j in tau_hat[-1])
{
  CI_p<-0
  CI_j<-NULL
  min_a<-1-j
  max_b<-length(Delta_hat)-j
  a<-0
  b<-0
  P_value_i_upper[j+a:b]<-(abs((1-CI_p)+P_value_i_upper[j+a:b])+abs((1-CI_p)-P_value_i_upper[j+a:b]))/2
  CI_p<-sum(Delta_pp[j+a:b])
  CI_j<-cbind(CI_j,c(j+a,j+b,CI_p))
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
      CI_j<-cbind(CI_j,c(j+a,j+b,CI_p))
      P_value_i_lower[j+a:b]<-(abs((1-CI_p)+P_value_i_lower[j+a:b])+abs((1-CI_p)-P_value_i_lower[j+a:b]))/2
    }else
    {
      CI_p<-1
      CI_j<-cbind(CI_j,c(j+a,j+b,CI_p))
    }
  }
  CI_tau[[which(tau_hat==j)-1]]<-CI_j
}


png("Delta_Summary3.png",600,480)
par(mar=c(4,6,1,1),mfrow=c(1,1))
plot(2:T_max,Delta_mean[-1],xlim=c(0,T_max),ylim = c(0,1.3),main = "",ylab="",xlab = "",type="b",pch=18,lty=1,col=2,axes = F, xaxs="i")
#mtext("Change point detection",side=2,cex=2)
mtext("Posterior Probability",side=2,line = 3,cex=2,at=0.5)
rect(0,-0.02,T_max,1,lwd=1.5)
axis(2,at=0:5/5,labels = 0:5/5,cex.axis=1.5)
for(i in (which(Delta_final==1))[-1])
{
  segments(i,0,i,1.3,col=1,lty=2)
  mtext(format(as.Date(as.character(X$Date[i]),format = "%d/%m/%Y"),"%b"),side = 1,at=i,line=0,cex=2)
  mtext(format(as.Date(as.character(X$Date[i]),format = "%d/%m/%Y"),"%d"),side = 1,at=i,line=2,cex=2)
}
alpha_list<-c(0.9,0.95,0.99)
for(alpha in 1:3)
{
  abline(h=1+0.1*alpha,lwd=0.8,lty=3)
  mtext(paste(c(as.character(alpha_list[alpha]*100),"% HPD"),collapse = ""),side = 2,at=1+0.1*alpha,las=2,line = 0,cex=1.5)
  for(aa in 1:length(CI_tau))
  {
    hai<-min(which(CI_tau[[aa]][3,]>=alpha_list[alpha]))
    segments(CI_tau[[aa]][1,hai],1+0.1*alpha,CI_tau[[aa]][2,hai],1+0.1*alpha,lwd=2,col="blue")
    segments(CI_tau[[aa]][1,hai],1+0.1*alpha-0.02,CI_tau[[aa]][1,hai],1+0.1*alpha+0.02,lwd=2,col="blue")
    segments(CI_tau[[aa]][2,hai],1+0.1*alpha-0.02,CI_tau[[aa]][2,hai],1+0.1*alpha+0.02,lwd=2,col="blue")
  }
  
}
dev.off()
b1<-colMeans(1/B_all)#,1/colMeans(B_all),
r2<-colMeans(R_all/(1+R_all))

png("beta_Summary3.png",600,600)
par(mar=c(6,6,2,1),mfrow=c(1,1))
plot(1:T_max,(b1),xlim=c(0,T_max),main = "",ylab="",xlab = "",type="b",pch=18,lty=1,col=2,xaxt="n",yaxt="n",cex=2,lwd=2)
mtext("Transmission rate",side=2,line = 3,cex=3)
axis(2,at=1:6/20,labels = 1:6/20,cex.axis=1.5)
for(i in (which(Delta_final==1))[-1])
{
  segments(i,0,i,1.3,col=1,lty=2)
  mtext(format(as.Date(as.character(X$Date[i]),format = "%d/%m/%Y"),"%b"),side = 1,at=i,line=1,cex=2.2)
  mtext(format(as.Date(as.character(X$Date[i]),format = "%d/%m/%Y"),"%d"),side = 1,at=i,line=3,cex=2.2)
}
dev.off()
png("gamma_Summary3.png",600,600)
par(mar=c(6,6,2,1),mfrow=c(1,1))
plot(1:T_max,(r2),xlim=c(0,T_max),main = "",ylab="",xlab = "",type="b",pch=18,lty=1,col="blue",xaxt="n",yaxt="n",cex=2,lwd=2)
mtext("Removal rate",side=2,line = 3,cex=3)
axis(2,at=1:4/10,labels = 1:4/10,cex.axis=1.5)
for(i in (which(Delta_final==1))[-1])
{
  segments(i,0,i,1.3,col=1,lty=2)
  mtext(format(as.Date(as.character(X$Date[i]),format = "%d/%m/%Y"),"%b"),side = 1,at=i,line=1,cex=2.2)
  mtext(format(as.Date(as.character(X$Date[i]),format = "%d/%m/%Y"),"%d"),side = 1,at=i,line=3,cex=2.2)
}
dev.off()