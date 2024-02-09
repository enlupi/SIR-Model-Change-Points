papapa<-function(jishu)
{
  Result<-NULL
  for(iter in 1:jishu)
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
    
    Burnin<-5000
    AAA<-1000
    thining<-10
    
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
      print(iter)
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
        print((iter-1)*thining+ai+Burnin)
      }
      Beta_all<-rbind(Beta_all,beta_hat)
      B_all<-rbind(B_all,b_hat[Stage_hat])
      
      Gamma_all<-rbind(Gamma_all,gamma_hat)
      R_all<-rbind(R_all,r_hat[Stage_hat])
      
      Delta_all<-rbind(Delta_all,Delta_hat)
    }
    
    
    
    
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
    Result<-rbind(Result,Delta_final,Delta_mean,colMeans(1/B_all),1/colMeans(B_all),colMeans(R_all/(1+R_all)),colMeans(R_all)/(1+colMeans(R_all)))
  }
  return(Result)
}

library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

wosai<- foreach(jishu=rep(1,10),.combine='rbind') %dopar% papapa(jishu)
stopCluster(cl)
write.csv(wosai,"Day100IIb_a.csv")

library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

wosai<- foreach(jishu=rep(1,10),.combine='rbind') %dopar% papapa(jishu)
stopCluster(cl)
write.csv(wosai,"Day100IIb_b.csv")

library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

wosai<- foreach(jishu=rep(1,10),.combine='rbind') %dopar% papapa(jishu)
stopCluster(cl)
write.csv(wosai,"Day100IIb_c.csv")

library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

wosai<- foreach(jishu=rep(1,10),.combine='rbind') %dopar% papapa(jishu)
stopCluster(cl)
write.csv(wosai,"Day100IIb_d.csv")

library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

wosai<- foreach(jishu=rep(1,10),.combine='rbind') %dopar% papapa(jishu)
stopCluster(cl)
write.csv(wosai,"Day100IIb_e.csv")

library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

wosai<- foreach(jishu=rep(1,10),.combine='rbind') %dopar% papapa(jishu)
stopCluster(cl)
write.csv(wosai,"Day100IIb_f.csv")

library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

wosai<- foreach(jishu=rep(1,10),.combine='rbind') %dopar% papapa(jishu)
stopCluster(cl)
write.csv(wosai,"Day100IIb_g.csv")

library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

wosai<- foreach(jishu=rep(1,10),.combine='rbind') %dopar% papapa(jishu)
stopCluster(cl)
write.csv(wosai,"Day100IIb_h.csv")

library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

wosai<- foreach(jishu=rep(1,10),.combine='rbind') %dopar% papapa(jishu)
stopCluster(cl)
write.csv(wosai,"Day100IIb_i.csv")

library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

wosai<- foreach(jishu=rep(1,10),.combine='rbind') %dopar% papapa(jishu)
stopCluster(cl)
write.csv(wosai,"Day100IIb_j.csv")

