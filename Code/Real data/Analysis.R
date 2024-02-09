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

########################## Hyperparameter
p_a<-1/T_max^3
p_b<-2-p_a

b_shape<-0.1
gamma_b_shape<-gamma(b_shape)
b_rate<-0.1

r_shape<-0.1
gamma_r_shape<-gamma(r_shape)
r_rate<-0.1
########################## Initialization
Delta_hat<-c(1,rep(0,T_max-2),0)
index_change<-1
Stage_hat<-cumsum(Delta_hat)
K_hat<-sum(Delta_hat)

b_hat<-rep(0,K_hat)
r_hat<-rep(0,K_hat)
for(i in 1:K_hat)
{
  b_hat[i]<-rgamma(1,shape=b_shape,rate=b_rate)
  r_hat[i]<-rgamma(1,shape=r_shape,rate=r_rate)
}

beta_hat<-rep(0,T_max)
gamma_hat<-rep(0,T_max)
contact_hat<-rep(0,T_max)
for(t in 1:T_max)
{
  beta_hat[t]<-rexp(1,rate = b_hat[Stage_hat[t]])
  gamma_hat[t]<-rbeta(1,shape1 = r_hat[Stage_hat[t]],shape2 = 1)
  
  lambda_t<-beta_hat[t]*PI_obs[t]
  p_upper<-1-ppois(0,lambda_t)
  pp<-runif(N_infect_obs[t])*p_upper+(1-p_upper)
  contact_hat[t]<-sum(qpois(pp,lambda_t))
}
for(t in 1:T_max)
{
  beta_hat[t]<-rgamma(1,shape = contact_hat[t]+1,rate = b_hat[Stage_hat[t]]+PI_obs[t]*S_obs[t])
  gamma_hat[t]<-rbeta(1,shape1 = N_recovery_obs[t]+r_hat[Stage_hat[t]],shape2 = 1+I_obs[t]-N_recovery_obs[t])
  
  lambda_t<-beta_hat[t]*PI_obs[t]
  p_upper<-1-ppois(0,lambda_t)
  pp<-runif(N_infect_obs[t])*p_upper+(1-p_upper)
  contact_hat[t]<-sum(qpois(pp,lambda_t))
}

max_k<-1+(T_max-1)

Accept_or_not<-NULL
for(iter in 1:Burnin)
{
  change_type<-sample((-1):1,1)
  # change_type<-1
  
  if(K_hat==1)
  {
    change_type<-1
  }
  
  if(K_hat==max_k)
  {
    change_type<--1
  }
  
  if(change_type==0)
  {
    change_position_all<-which(abs(Delta_hat[2:(T_max-1)]-Delta_hat[1+2:(T_max-1)])==1)+1
    change_position<-change_position_all[sample(length(change_position_all),1)]
    
    Delta_hat_candidate<-Delta_hat
    Delta_hat_candidate[change_position+0:1]<-Delta_hat_candidate[change_position+1:0]
    Stage_hat_candidate<-cumsum(Delta_hat_candidate)
    index_change_candidate<-which(Delta_hat_candidate==1)
    change_position_all_candidate<-which(abs(Delta_hat_candidate[2:(T_max-1)]-Delta_hat_candidate[1+2:(T_max-1)])==1)+1
    
    phase<-c(Stage_hat[change_position],Stage_hat_candidate[change_position])
    logp_original<-0
    logp_candidate<-0
    for(i in phase)
    {
      L_i_original<-which(Stage_hat==i)
      logp_original<-logp_original+sum(log(1:length(L_i_original)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_original))*log(b_rate+sum(beta_hat[L_i_original]))
      logp_original<-logp_original+sum(log(1:length(L_i_original)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_original))*log(r_rate+sum(-log(gamma_hat[L_i_original])))
      logp_original<-logp_original-sum(log(1:length(L_i_original)))
      
      L_i_candidate<-which(Stage_hat_candidate==i)
      logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_candidate))*log(b_rate+sum(beta_hat[L_i_candidate]))
      logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_candidate))*log(r_rate+sum(-log(gamma_hat[L_i_candidate])))
      logp_candidate<-logp_candidate-sum(log(1:length(L_i_candidate)))
    }
    
    ratio<-exp(min(c(0,logp_candidate-logp_original+log(length(change_position_all)/length(change_position_all_candidate)))))
    cxx<-rbinom(1,1,ratio)
    if(cxx==1)
    {
      Delta_hat<-Delta_hat_candidate
      Stage_hat<-Stage_hat_candidate
      index_change<-which(Delta_hat==1)
    }
  }else
  {
    if(change_type==1)
    {
      change_position_all<-which(Delta_hat[2:(T_max)]==0)+1
      change_position<-change_position_all[sample(length(change_position_all),1)]
      
      Delta_hat_candidate<-Delta_hat
      Delta_hat_candidate[change_position]<-1-Delta_hat_candidate[change_position]
      Stage_hat_candidate<-cumsum(Delta_hat_candidate)
      index_change_candidate<-which(Delta_hat_candidate==1)
      change_position_all_candidate<-which(Delta_hat_candidate[2:(T_max)]==1)+1
      
      phase_original<-Stage_hat[change_position]
      phase_candidate<-c(Stage_hat[change_position],Stage_hat_candidate[change_position])
      
      logp_original<-0
      logp_candidate<-0
      for(i in phase_original)
      {
        L_i_original<-which(Stage_hat==i)
        logp_original<-logp_original+sum(log(1:length(L_i_original)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_original))*log(b_rate+sum(beta_hat[L_i_original]))
        logp_original<-logp_original+sum(log(1:length(L_i_original)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_original))*log(r_rate+sum(-log(gamma_hat[L_i_original])))
        logp_original<-logp_original-sum(log(1:length(L_i_original)))
      }
      
      for(i in phase_candidate)
      {
        L_i_candidate<-which(Stage_hat_candidate==i)
        logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_candidate))*log(b_rate+sum(beta_hat[L_i_candidate]))
        logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_candidate))*log(r_rate+sum(-log(gamma_hat[L_i_candidate])))
        logp_candidate<-logp_candidate-sum(log(1:length(L_i_candidate)))
      }
      logp_candidate<-logp_candidate+log(p_a/p_b)+log((3-2*(K_hat==1))*length(change_position_all)/(3-2*((K_hat+1)==max_k))/length(change_position_all_candidate))
      
      ratio<-exp(min(c(0,(logp_candidate-logp_original))))
      cxx<-rbinom(1,1,ratio)
      
      if(cxx==1)
      {
        Delta_hat<-Delta_hat_candidate
        Stage_hat<-Stage_hat_candidate
      }
    }else
    {
      # Delta_hat_candidate<-Delta_hat
      change_position_all<-which(Delta_hat[2:(T_max)]==1)+1
      change_position<-change_position_all[sample(length(change_position_all),1)]
      
      Delta_hat_candidate<-Delta_hat
      Delta_hat_candidate[change_position]<-1-Delta_hat_candidate[change_position]
      Stage_hat_candidate<-cumsum(Delta_hat_candidate)
      index_change_candidate<-which(Delta_hat_candidate==1)
      change_position_all_candidate<-which(Delta_hat_candidate[2:(T_max)]==0)+1
      
      phase_original<-c(Stage_hat[change_position],Stage_hat_candidate[change_position])
      phase_candidate<-c(Stage_hat_candidate[change_position])
      
      logp_original<-0
      logp_candidate<-0
      for(i in phase_original)
      {
        L_i_original<-which(Stage_hat==i)
        logp_original<-logp_original+sum(log(1:length(L_i_original)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_original))*log(b_rate+sum(beta_hat[L_i_original]))
        logp_original<-logp_original+sum(log(1:length(L_i_original)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_original))*log(r_rate+sum(-log(gamma_hat[L_i_original])))
        logp_original<-logp_original-sum(log(1:length(L_i_original)))
      }
      
      for(i in phase_candidate)
      {
        L_i_candidate<-which(Stage_hat_candidate==i)
        logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_candidate))*log(b_rate+sum(beta_hat[L_i_candidate]))
        logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_candidate))*log(r_rate+sum(-log(gamma_hat[L_i_candidate])))
        logp_candidate<-logp_candidate-sum(log(1:length(L_i_candidate)))
      }
      logp_candidate<-logp_candidate-log(p_a/p_b)+log((3-2*(K_hat==max_k))*length(change_position_all)/(3-2*((K_hat-1)==1))/length(change_position_all_candidate))
      
      ratio<-exp(min(c(0,(logp_candidate-logp_original))))
      cxx<-rbinom(1,1,ratio)
      
      if(cxx==1)
      {
        Delta_hat<-Delta_hat_candidate
        Stage_hat<-Stage_hat_candidate
      }
    }
  }
  
  K_hat<-sum(Delta_hat)
  
  b_hat<-rep(0,K_hat)
  r_hat<-rep(0,K_hat)
  for(i in 1:K_hat)
  {
    L_i<-which(Stage_hat==i)
    b_hat[i]<-rgamma(1,shape=b_shape+(length(L_i)),rate=b_rate+sum(beta_hat[L_i]))
    r_hat[i]<-rgamma(1,shape=r_shape+(length(L_i)),rate=r_rate+sum(-log(gamma_hat[L_i])))
  }
  for(t in 1:T_max)
  {
    beta_hat[t]<-rgamma(1,shape = contact_hat[t]+1,rate = b_hat[Stage_hat[t]]+PI_obs[t]*S_obs[t])
    gamma_hat[t]<-rbeta(1,shape1 = N_recovery_obs[t]+r_hat[Stage_hat[t]],shape2 = 1+I_obs[t]-N_recovery_obs[t])
    
    lambda_t<-beta_hat[t]*PI_obs[t]
    p_upper<-1-ppois(0,lambda_t)
    pp<-runif(N_infect_obs[t])*p_upper+(1-p_upper)
    contact_hat[t]<-sum(qpois(pp,lambda_t))
  }
  # print(iter)
}
Beta_all<-NULL
B_all<-NULL

Gamma_all<-NULL
R_all<-NULL

Delta_all<-NULL

for(iter in 1:AAA)
{
  for(ai in 1:thining)
  {
    change_type<-sample((-1):1,1)
    # change_type<-1
    
    if(K_hat==1)
    {
      change_type<-1
    }
    
    if(K_hat==max_k)
    {
      change_type<--1
    }
    
    if(change_type==0)
    {
      change_position_all<-which(abs(Delta_hat[2:(T_max-1)]-Delta_hat[1+2:(T_max-1)])==1)+1
      change_position<-change_position_all[sample(length(change_position_all),1)]
      
      Delta_hat_candidate<-Delta_hat
      Delta_hat_candidate[change_position+0:1]<-Delta_hat_candidate[change_position+1:0]
      Stage_hat_candidate<-cumsum(Delta_hat_candidate)
      index_change_candidate<-which(Delta_hat_candidate==1)
      change_position_all_candidate<-which(abs(Delta_hat_candidate[2:(T_max-1)]-Delta_hat_candidate[1+2:(T_max-1)])==1)+1
      
      phase<-c(Stage_hat[change_position],Stage_hat_candidate[change_position])
      logp_original<-0
      logp_candidate<-0
      for(i in phase)
      {
        L_i_original<-which(Stage_hat==i)
        logp_original<-logp_original+sum(log(1:length(L_i_original)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_original))*log(b_rate+sum(beta_hat[L_i_original]))
        logp_original<-logp_original+sum(log(1:length(L_i_original)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_original))*log(r_rate+sum(-log(gamma_hat[L_i_original])))
        logp_original<-logp_original-sum(log(1:length(L_i_original)))
        
        L_i_candidate<-which(Stage_hat_candidate==i)
        logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_candidate))*log(b_rate+sum(beta_hat[L_i_candidate]))
        logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_candidate))*log(r_rate+sum(-log(gamma_hat[L_i_candidate])))
        logp_candidate<-logp_candidate-sum(log(1:length(L_i_candidate)))
      }
      
      ratio<-exp(min(c(0,logp_candidate-logp_original+log(length(change_position_all)/length(change_position_all_candidate)))))
      cxx<-rbinom(1,1,ratio)
      if(cxx==1)
      {
        Delta_hat<-Delta_hat_candidate
        Stage_hat<-Stage_hat_candidate
        index_change<-which(Delta_hat==1)
      }
    }else
    {
      if(change_type==1)
      {
        change_position_all<-which(Delta_hat[2:(T_max)]==0)+1
        change_position<-change_position_all[sample(length(change_position_all),1)]
        
        Delta_hat_candidate<-Delta_hat
        Delta_hat_candidate[change_position]<-1-Delta_hat_candidate[change_position]
        Stage_hat_candidate<-cumsum(Delta_hat_candidate)
        index_change_candidate<-which(Delta_hat_candidate==1)
        change_position_all_candidate<-which(Delta_hat_candidate[2:(T_max)]==1)+1
        
        phase_original<-Stage_hat[change_position]
        phase_candidate<-c(Stage_hat[change_position],Stage_hat_candidate[change_position])
        
        logp_original<-0
        logp_candidate<-0
        for(i in phase_original)
        {
          L_i_original<-which(Stage_hat==i)
          logp_original<-logp_original+sum(log(1:length(L_i_original)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_original))*log(b_rate+sum(beta_hat[L_i_original]))
          logp_original<-logp_original+sum(log(1:length(L_i_original)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_original))*log(r_rate+sum(-log(gamma_hat[L_i_original])))
          logp_original<-logp_original-sum(log(1:length(L_i_original)))
        }
        
        for(i in phase_candidate)
        {
          L_i_candidate<-which(Stage_hat_candidate==i)
          logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_candidate))*log(b_rate+sum(beta_hat[L_i_candidate]))
          logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_candidate))*log(r_rate+sum(-log(gamma_hat[L_i_candidate])))
          logp_candidate<-logp_candidate-sum(log(1:length(L_i_candidate)))
        }
        logp_candidate<-logp_candidate+log(p_a/p_b)+log((3-2*(K_hat==1))*length(change_position_all)/(3-2*((K_hat+1)==max_k))/length(change_position_all_candidate))
        
        ratio<-exp(min(c(0,(logp_candidate-logp_original))))
        cxx<-rbinom(1,1,ratio)
        
        if(cxx==1)
        {
          Delta_hat<-Delta_hat_candidate
          Stage_hat<-Stage_hat_candidate
        }
      }else
      {
        # Delta_hat_candidate<-Delta_hat
        change_position_all<-which(Delta_hat[2:(T_max)]==1)+1
        change_position<-change_position_all[sample(length(change_position_all),1)]
        
        Delta_hat_candidate<-Delta_hat
        Delta_hat_candidate[change_position]<-1-Delta_hat_candidate[change_position]
        Stage_hat_candidate<-cumsum(Delta_hat_candidate)
        index_change_candidate<-which(Delta_hat_candidate==1)
        change_position_all_candidate<-which(Delta_hat_candidate[2:(T_max)]==0)+1
        
        phase_original<-c(Stage_hat[change_position],Stage_hat_candidate[change_position])
        phase_candidate<-c(Stage_hat_candidate[change_position])
        
        logp_original<-0
        logp_candidate<-0
        for(i in phase_original)
        {
          L_i_original<-which(Stage_hat==i)
          logp_original<-logp_original+sum(log(1:length(L_i_original)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_original))*log(b_rate+sum(beta_hat[L_i_original]))
          logp_original<-logp_original+sum(log(1:length(L_i_original)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_original))*log(r_rate+sum(-log(gamma_hat[L_i_original])))
          logp_original<-logp_original-sum(log(1:length(L_i_original)))
        }
        
        for(i in phase_candidate)
        {
          L_i_candidate<-which(Stage_hat_candidate==i)
          logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+b_shape))-log(gamma_b_shape)+b_shape*log(b_rate)-(b_shape+length(L_i_candidate))*log(b_rate+sum(beta_hat[L_i_candidate]))
          logp_candidate<-logp_candidate+sum(log(1:length(L_i_candidate)-1+r_shape))-log(gamma_r_shape)+r_shape*log(r_rate)-(r_shape+length(L_i_candidate))*log(r_rate+sum(-log(gamma_hat[L_i_candidate])))
          logp_candidate<-logp_candidate-sum(log(1:length(L_i_candidate)))
        }
        logp_candidate<-logp_candidate-log(p_a/p_b)+log((3-2*(K_hat==max_k))*length(change_position_all)/(3-2*((K_hat-1)==1))/length(change_position_all_candidate))
        
        ratio<-exp(min(c(0,(logp_candidate-logp_original))))
        cxx<-rbinom(1,1,ratio)
        
        if(cxx==1)
        {
          Delta_hat<-Delta_hat_candidate
          Stage_hat<-Stage_hat_candidate
        }
      }
    }
    
    K_hat<-sum(Delta_hat)
    
    b_hat<-rep(0,K_hat)
    r_hat<-rep(0,K_hat)
    for(i in 1:K_hat)
    {
      L_i<-which(Stage_hat==i)
      b_hat[i]<-rgamma(1,shape=b_shape+(length(L_i)),rate=b_rate+sum(beta_hat[L_i]))
      r_hat[i]<-rgamma(1,shape=r_shape+(length(L_i)),rate=r_rate+sum(-log(gamma_hat[L_i])))
    }
    for(t in 1:T_max)
    {
      beta_hat[t]<-rgamma(1,shape = contact_hat[t]+1,rate = b_hat[Stage_hat[t]]+PI_obs[t]*S_obs[t])
      gamma_hat[t]<-rbeta(1,shape1 = N_recovery_obs[t]+r_hat[Stage_hat[t]],shape2 = 1+I_obs[t]-N_recovery_obs[t])
      
      lambda_t<-beta_hat[t]*PI_obs[t]
      p_upper<-1-ppois(0,lambda_t)
      pp<-runif(N_infect_obs[t])*p_upper+(1-p_upper)
      contact_hat[t]<-sum(qpois(pp,lambda_t))
    }
    # print((iter-1)*thining+ai+Burnin)
  }
  Beta_all<-rbind(Beta_all,beta_hat)
  B_all<-rbind(B_all,b_hat[Stage_hat])
  
  Gamma_all<-rbind(Gamma_all,gamma_hat)
  R_all<-rbind(R_all,r_hat[Stage_hat])
  
  Delta_all<-rbind(Delta_all,Delta_hat)
}
write.csv(Accept_or_not,"Accept_or_not.csv",row.names = F)
write.csv(Beta_all,"Beta_all.csv",row.names = F)
write.csv(B_all,"B_all.csv",row.names = F)
write.csv(Gamma_all,"Gamma_all.csv",row.names = F)
write.csv(R_all,"R_all.csv",row.names = F)
write.csv(Delta_all,"Delta_all.csv",row.names = F)