import numpy as np
from scipy.special import gamma as gammaFunc
from scipy.stats import gamma, beta
import tqdm



# New delta proposal
def propose_delta(Delta_hat, K_hat, T_max=100, verbose=False):

    ########## change type
    ## add–delete–swap - proposal step
    # -1 delete
    # 0 swap
    # +1 add
    if K_hat==1:
        change_type = 1
        
    elif K_hat==T_max:
        change_type = -1

    else:
        change_type = np.random.choice([-1, 0, 1])
            
    Delta_hat_candidate = Delta_hat.copy()
    
    ########## proposal phase
    if change_type != 0:
        if change_type == 1:
            possible_change_indices = np.where(Delta_hat[1:] == 0)[0]+1
        if change_type == -1:
            possible_change_indices = np.where(Delta_hat[1:] == 1)[0]+1
            
        index_to_change = np.random.choice(possible_change_indices)
            
        Delta_hat_candidate[index_to_change] = 1-Delta_hat_candidate[index_to_change]
            
    else:
        possible_change_indices = np.where(np.abs(Delta_hat[1:-1] - Delta_hat[2:]) == 1)[0]+1
        index_to_change = np.random.choice(possible_change_indices)

        Delta_hat_candidate[index_to_change + np.array([0, 1])] = Delta_hat_candidate[index_to_change + np.array([1, 0])]

    Stage_hat_candidate = np.cumsum(Delta_hat_candidate, dtype=int)-1
    K_hat_candidate = np.sum(Delta_hat_candidate, dtype=int)
    
    if verbose:
        print("Change:", change_type)
        print("Original:",  np.where(Delta_hat==1)[0])
        print("Candidate:", np.where(Delta_hat_candidate==1)[0])
        print("---------------------------")

    return Delta_hat_candidate, Stage_hat_candidate, K_hat_candidate



# delta sampling
def sample_delta_1(Delta_hat, Stage_hat, K_hat, beta_hat, gamma_hat, 
                   p=0.01, b_shape=0.1, b_rate=0.1, r_shape=0.1, r_rate=0.1, T_max=100, verbose=False):
    
    # propose new delta
    Delta_hat_candidate, Stage_hat_candidate, K_hat_candidate = propose_delta(Delta_hat, K_hat, T_max, verbose)


    ########## Metropolis-Hastings Phase
    # Prior log-ratio: pi(d*)/pi(d(g))
    log_pi = (K_hat_candidate-K_hat)*np.log(p/(1-p))
        
    # Likelihood log-ratio: pi(beta(g), gamma(g)|d*)/pi(beta(g), gamma(g)|d(g))
    log1_candidate = 0
    log2_candidate = 0
        
    for k in range(K_hat_candidate):
        ind_k = np.where(Stage_hat_candidate == k)[0]
        log1_candidate += np.log(gammaFunc(b_shape+len(ind_k))) - (b_shape+len(ind_k))*np.log(b_rate+np.sum(beta_hat[ind_k]))
        log2_candidate += np.log(gammaFunc(r_shape+len(ind_k))) - (r_shape+len(ind_k))*np.log(r_rate+np.sum(-np.log(gamma_hat[ind_k])))
        
    log1_original = 0
    log2_original = 0
                                     
    for k in range(K_hat):
        ind_k = np.where(Stage_hat == k)[0]
        log1_original += np.log(gammaFunc(b_shape+len(ind_k))) - (b_shape+len(ind_k))*np.log(b_rate+np.sum(beta_hat[ind_k]))
        log2_original += np.log(gammaFunc(r_shape+len(ind_k))) - (r_shape+len(ind_k))*np.log(r_rate+np.sum(-np.log(gamma_hat[ind_k])))
        
    log_L = (log1_candidate+log2_candidate) - (log1_original+log2_original)     
    
    # Jump probability log-ratio: J(d(g)|d*)/J(d*|d(g))
    JJ = 0
    if K_hat == K_hat_candidate:
        JJ = 1 
    elif ([K_hat_candidate, K_hat] == [1, 2] or [K_hat_candidate, K_hat] == [T_max, T_max-1]):
        JJ = 3/(T_max-1)
    elif ([K_hat_candidate, K_hat] == [2, 1] or [K_hat_candidate, K_hat] == [T_max-1, T_max]):
        JJ = (T_max-1)/3
    elif (K_hat_candidate-K_hat) == -1 and K_hat_candidate != 1 and K_hat_candidate != (T_max-1):
        JJ = (K_hat-1)/(T_max - K_hat_candidate)
    elif K_hat_candidate - K_hat == 1 and K_hat_candidate != 2 and K_hat_candidate != T_max:
        JJ = (T_max-K_hat)/(K_hat_candidate-1)
    log_JJ = np.log(JJ)

    
    # Metropolis-Hastings Ratio
    log_mMH = log_L+log_pi+log_JJ
    ratio = np.exp(min(0, log_mMH))

    if verbose:
        print("pi:", log_pi, np.exp(log_pi))
        print("likelihood:", log_L, np.exp(log_L))
        print("JJ:", log_JJ, JJ)
        print("m_MH:", log_mMH, np.exp(log_mMH))
    
    cxx = np.random.binomial(1, ratio) 
    if cxx == 1:
        Delta_hat = Delta_hat_candidate
        Stage_hat = Stage_hat_candidate
        K_hat     = K_hat_candidate
        
    return Delta_hat, Stage_hat, K_hat



# Get MCMC chain using Gibbs sampling 
def gibbs_sampling(data, samples=1000, T_max=100, burnin=5000, thinning=10, 
                   delta_sample='paper', p=0.01, chg_pt=None, verbose=False):
    
    ############# data
    # Data is expected to be a dataframe with 101 rows (100 steps + initial one)
    # OBS: the initial row should display the starting values 
    #      for Susceptible (S), Infected (I), and Recovered (R), with random values in the 
    #      columns representing changes (deltas) over time.
    
    I       = data["I"].values
    S       = data["S"].values
    PI      = data["PI"].values
    Delta_I = data["Delta_I"].values[1:]
    Delta_R = data["Delta_R"].values[1:]
    
    
    #######################
    ##  HYPERPARAMETERS  ##
    #######################
    
    b_shape = 0.1
    b_rate  = 0.1
    r_shape = 0.1
    r_rate  = 0.1

    ######################
    ##  INITIALIZATION  ##
    ######################
    
    ##### delta
    Delta_hat = np.zeros(shape=T_max, dtype=int)

    if delta_sample=='fixed':
        Delta_hat[chg_pt] = 1
    else:
        Delta_hat[0] = 1
        
    Stage_hat = np.cumsum(Delta_hat, dtype=int)-1
    K_hat     = np.sum(Delta_hat, dtype=int)
                        
    
    ##### b and r
    b_hat    = gamma.rvs(a=b_shape, scale=1/b_rate, size=K_hat)
    r_hat    = gamma.rvs(a=r_shape, scale=1/r_rate, size=K_hat)
    

    ##### beta and gamma
    
    chi_hat   = beta.rvs(a=1+Delta_I,                b=b_hat[Stage_hat]/PI[:-1]+S[:-1]-Delta_I)
    beta_hat  = -np.log(1-chi_hat)/PI[:-1] 
    gamma_hat = beta.rvs(a=r_hat[Stage_hat]+Delta_R, b=1+I[:-1]-Delta_R)
    gamma_hat[gamma_hat == 0] = 10e-7
    
    print("Initialization Complete\n")
    if verbose:
        print("Delta_hat:",Delta_hat)
        print("b_hat:",b_hat)
        print("r_hat:",r_hat)
        print("beta_hat:",beta_hat)
        print("gamma_hat:",gamma_hat)
        print("----------------------------------------------------")
        print("----------------------------------------------------")
    
    
    ################
    ##  SAMPLING  ##
    ################
    
    Delta_all = []
    Stage_all = []
    b_all     = []
    r_all     = []
    beta_all  = []
    gamma_all = []
    
    for step in tqdm.tqdm(range(burnin+thinning*samples)):
        
        ##### delta sampling

        if delta_sample=='paper':
            Delta_hat, Stage_hat, K_hat = sample_delta_1(Delta_hat, Stage_hat, K_hat, beta_hat, gamma_hat, 
                                                         p, b_shape, b_rate, r_shape, r_rate, T_max)
        #elif delta_sample=='fixed':
        #   no change 
        
        ##### b and r sampling
        b_hat = np.zeros(K_hat)
        r_hat = np.zeros(K_hat)
        for i in range(K_hat):
            L_i = np.where(Stage_hat == i)[0]
            b_hat[i] = gamma.rvs(a=(b_shape + len(L_i)), scale=1/(b_rate + np.sum(beta_hat[L_i])))
            r_hat[i] = gamma.rvs(a=(r_shape + len(L_i)), scale=1/(r_rate + np.sum(-np.log(gamma_hat[L_i]))))
        
        ##### beta and gamma sampling
        chi_hat   = beta.rvs(a=1+Delta_I,b=b_hat[Stage_hat]/PI[:-1]+S[:-1]-Delta_I)
        beta_hat  = -np.log(1-chi_hat)/PI[:-1] 
        
        gamma_hat = beta.rvs(scale=1,a=r_hat[Stage_hat]+Delta_R, b=1+I[:-1]-Delta_R)
        gamma_hat[gamma_hat == 0] = 10e-7
        
        if verbose:
            if step % 100 == 0 and step != 0:
                print("\nStep:",step)
                print("Delta_hat:",Delta_hat)
                print("b_hat:",b_hat)
                print("r_hat:",r_hat)
                print("beta_hat:",beta_hat)
                print("gamma_hat:",gamma_hat)
                print("----------------------------------------------------")
            
        Delta_all.append(Delta_hat)
        Stage_all.append(Stage_hat)
        b_all.append(b_hat)
        r_all.append(r_hat)
        beta_all.append(beta_hat)
        gamma_all.append(gamma_hat)
  
    # DataFrame generation
    MCMC_chain = {
        'Delta': Delta_all[burnin::thinning],
        'Stage': Stage_all[burnin::thinning],
        'b': b_all[burnin::thinning],
        'r': r_all[burnin::thinning],
        'beta': beta_all[burnin::thinning],
        'gamma': gamma_all[burnin::thinning]
        }
    
    return MCMC_chain