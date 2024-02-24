import numpy as np
import itertools



# Loss function for delta estimation
def compute_loss(delta, Q_tt):
    stage = np.cumsum(delta, dtype=int)-1
    Q_candidate = (stage[:, np.newaxis] == stage)

    return np.sum(np.abs(Q_candidate - Q_tt))



# Estimators for delta, beta and gamma
def get_estimators(chain, T_max=100, samples=1000):

    ##### Delta
    # get Bayes estimator
    stage = np.array(chain["Stage"])
    masks = (stage[:, :, np.newaxis] == stage[:, np.newaxis, :])
    Q_ttprime = np.mean(masks, axis=0)

    Delta_final = np.zeros(shape=T_max, dtype=int)
    Delta_final[0] = 1
    current_loss = compute_loss(Delta_final, Q_ttprime)

    check_add_drop = True
    check_swap = True
    while(check_add_drop==True or check_swap==True):
        # check add or drop
        candidates_loss = np.zeros(shape=(T_max))
        for i in range(1, T_max):
            Delta_candidate = Delta_final.copy()
            Delta_candidate[i] = 1-Delta_candidate[i]
            candidates_loss[i] = compute_loss(Delta_candidate, Q_ttprime)
        index_min = np.argmin(candidates_loss[1:])+1
        if candidates_loss[index_min] < current_loss:
            current_loss = candidates_loss[index_min]
            Delta_final[index_min] = 1-Delta_final[index_min]
            check_add_drop = True
        else:
            check_add_drop = False

        # check swap
        if np.sum(Delta_final[1:]) in np.arange(1,T_max-1):
            possible_change_indices = np.where(np.abs(Delta_final[1:-1] - Delta_final[2:]) == 1)[0]+1
            candidates_loss = np.zeros(shape=(len(possible_change_indices)))
            for i, idx in enumerate(possible_change_indices):
                Delta_candidate = Delta_final.copy()
                Delta_candidate[idx+np.array([0,1])] = Delta_candidate[idx+np.array([1,0])]
                candidates_loss[i] = compute_loss(Delta_candidate, Q_ttprime)
            index_min = np.argmin(candidates_loss)
            if candidates_loss[index_min] < current_loss:
                current_loss = candidates_loss[index_min]
                Delta_candidate[possible_change_indices[index_min]+np.array([0,1])] = Delta_candidate[possible_change_indices[index_min]+np.array([1,0])]
                check_swap = True
            else:
                check_swap = False

    ##### beta and gamma
                
    b = np.array([[chain["b"][s][stage[s][t]] for t in range(T_max)] for s in range(samples)])
    r = np.array([[chain["r"][s][stage[s][t]] for t in range(T_max)] for s in range(samples)])
    beta_final  = 1/np.mean(b, axis=0)
    gamma_final = np.mean(r/(1+r), axis=0)

    return Delta_final, beta_final, gamma_final 



# Highest Posterior Density intervals for delta estimator
def get_HPD(delta, probs, alpha, T_max):
    idxs = np.where(delta[1:]==1)[0]+1
    HPD = np.zeros(shape=(len(idxs),2))
    for n, i in enumerate(idxs):
        l = 0 # length of the HPD interval candidate
        P = 0 # probability of the HPD interval candidate
        while True:
            i_left  = i-l
            i_right = i+l
            k_max = l

            if i_right > T_max-1:
                i_right = T_max-1
                k_max = T_max-1-i

            if i_left < 0:
                i_left = 0
                k_max = i

            p = np.array([np.sum(probs[i_left+k:i_left+k+l+1]) for k in range(k_max+1)]) # probabilities of each interval
                                                                                         # of length l containing i
            P = np.max(p)
            if P >= 1-alpha:
                k_star = np.argmax(p)
                HPD[n,0] = i_left+k_star
                HPD[n,1] = i_left+k_star+l
                break
            else:
                l += 1

    return HPD



# Compute change points probabilities and HPD intervals
def compute_probs(chain, delta_est, alpha, T_max=100):
    prob_delta = np.mean(np.array(chain["Delta"]), axis=0)
    HPD = get_HPD(delta_est, prob_delta, alpha, T_max)

    return prob_delta, HPD



# compute mean, mode, median
def comp_stat(chain, var, bins):
    array = np.array(chain[var])
    mean   = np.mean(array, axis=0)
    median = np.median(array, axis=0)
    mode = np.zeros(shape=len(mean))
    for i in range(len(mode)):
        x, bin = np.histogram(array[:,i], bins=bins)
        mode[i] = bin[np.argmax(x)] + (bin[1]-bin[0])/2

    return mean, mode, median



# Adjusted Rand Index
def comp_ARI(true, estim, T_max=100):

    comb = np.array(list(itertools.combinations(range(T_max), 2))).T
    true_mask = (true[comb[0]] == true[comb[1]]).astype(int)
    estim_mask = (estim[comb[0]] == estim[comb[1]]).astype(int)

    TP = np.mean(true_mask*estim_mask)
    FP = np.mean((1-true_mask)*estim_mask)
    FN = np.mean(true_mask*(1-estim_mask))
    TN = np.mean((1-true_mask)*(1-estim_mask))
    num = TP+TN-(TP+FP)*(TP+FN)-(TN+FP)*(TN+FN)
    den = 1-(TP+FP)*(TP+FN)-(TN+FP)*(TN+FN)

    ARI = num/den

    return ARI



# Mutual Information
def comp_MI(true, estim, T_max=100):

    _, counts = np.unique(true, return_counts=True)
    max_ = 0

    for n in counts:
        max_ += n/T_max*np.log(T_max/n)

    n_kkprime = np.histogram2d(true, estim, bins=(np.max(true)+1, np.max(estim)+1))[0]
    n_k = np.sum(n_kkprime, axis=1)
    n_kprime = np.sum(n_kkprime, axis=0)

    MI = np.sum(n_kkprime/T_max*np.log((n_kkprime+(n_kkprime==0))*T_max/np.outer(n_k, n_kprime)))

    return MI, max_