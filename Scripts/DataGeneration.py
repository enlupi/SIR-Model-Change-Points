import numpy as np
from scipy.stats import binom
import pandas as pd



# simulate dataset from initial conditions 
# can handle multiple datasets creations and scenarios at the same time
# returned arrays are in the shape: [n. of datasets; n. of scenarios; simulation time] 
def sim_dataset(chg_pt, scenarios, T_max, S0, I0, R0, n_datasets):

    N = S0 + I0 + R0
    n_sc = len(scenarios)

    # create array of transmission and removal rate parameters at each time step
    beta  = np.array([[scenarios[i][(chg_pt[i] <= t+1).sum()][0] for t in range(T_max)] for i in range(n_sc)])
    gamma = np.array([[scenarios[i][(chg_pt[i] <= t+1).sum()][1] for t in range(T_max)] for i in range(n_sc)])


    Delta_I = np.zeros(shape=(n_datasets, n_sc, T_max+1), dtype=np.int32)
    Delta_R = np.zeros(shape=(n_datasets, n_sc, T_max+1), dtype=np.int32)
    S       = np.zeros(shape=(n_datasets, n_sc, T_max+1), dtype=np.int32)
    I       = np.zeros(shape=(n_datasets, n_sc, T_max+1), dtype=np.int32)
    R       = np.zeros(shape=(n_datasets, n_sc, T_max+1), dtype=np.int32)

    Delta_I[:,:,0] = -1
    Delta_R[:,:,0] = -1
    S[:,:,0]       = S0
    I[:,:,0]       = I0
    R[:,:,0]       = R0

    for t in range(1, T_max+1):
        Delta_I[:,:,t] = binom.rvs(S[:,:,t-1], 1-np.exp(-beta[:,t-1]*I[:,:,t-1]/N))
        Delta_R[:,:,t] = binom.rvs(I[:,:,t-1], gamma[:,t-1])
        S[:,:,t]       = S[:,:,t-1] - Delta_I[:,:,t]
        I[:,:,t]       = I[:,:,t-1] + Delta_I[:,:,t] - Delta_R[:,:,t]
        R[:,:,t]       = R[:,:,t-1] + Delta_R[:,:,t]

    return Delta_I, Delta_R, S, I, R



# creates dataframe for particular scenario and dataset, cutting
# it short if the disease ends before the maximum simulation time
def sanitise_data(S, I, R, Delta_I, Delta_R, sc=0, d=0, N=1_000_000):

    PI = I.astype(np.float64)/N
    stop = S.shape[-1]-1
    zeros = np.where(I[d,sc]==0)[0]

    if len(zeros) > 0:
        stop = zeros[0]

    data = pd.DataFrame({
        'S': S[d,sc,:stop+1],
        'I': I[d,sc,:stop+1],
        'R': R[d,sc,:stop+1],
        'PI': PI[d,sc,:stop+1],
        'Delta_I': Delta_I[d,sc,:stop+1],
        'Delta_R': Delta_R[d,sc,:stop+1]
    })

    return data, stop