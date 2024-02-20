import numpy as np
from scipy.stats import binom, gamma, beta, expon, poisson, uniform, bernoulli
import scipy.stats as stats
from joblib import Parallel, delayed
from scipy.special import gamma as gammaFunc
import random

import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

import tqdm

import pickle 

import math

import itertools


def plot_SRI(S, I, R, sc=0, d=None, tot=1_000_000, time=100):

    if d is None:
        S = np.expand_dims(np.mean(S, axis=0), 0)
        I = np.expand_dims(np.mean(I, axis=0), 0)
        R = np.expand_dims(np.mean(R, axis=0), 0)
        d = 0

    S = S[d,sc]
    I = I[d,sc]
    R = R[d,sc]

    # plot
    fig, ax = plt.subplots()

    y = np.vstack([S, I, R])
    ax.stackplot(np.arange(time+1), y/tot, labels=["S","I","R"], alpha=0.8)

    ax.set_xlabel("Day")
    ax.set_ylabel("Proportion")
    ax.set_xticks(np.concatenate([[0], np.arange(25, time+1, 25)]))
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    plt.legend(loc="upper right")

    plt.show()
    
def plot_delta(prob_delta, delta_true, delta_est, HPD, T_max=100):

    fig, ax = plt.subplots(figsize=(10,5))
    ax.set_xticks(np.arange(5,T_max+1,10))
    ax.set_xlim(0, T_max+1)
    ymax = 1.3*np.max(prob_delta[1:])
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Day")
    ax.set_ylabel("Posterior probability", rotation=90)
    x = np.arange(2, T_max+1)

    ax.plot(x, prob_delta[1:], marker='.', markersize=4, color='red', linewidth=0.8)

    for i, x in enumerate(np.where(delta_est[1:] == 1)[0]+1):
        ax.vlines(x+1, ymin=0, ymax=ymax, color='red', linestyle='dashed', alpha=0.7)
        ax.axvspan(HPD[i,0]+1, HPD[i,1]+1, alpha=0.1, color='red')

    for x in np.where(delta_true[1:] == 1)[0]+1:
        ax.vlines(x+1, ymin=0, ymax=ymax, color='black', linestyle='dashed', alpha=0.7)
        
    l1, = plt.plot([-1], [-1], color='red', marker='.', markersize=4, label='Probability')
    l2, = plt.plot([-1], [-1], color='red', linestyle='dashed', alpha=0.7, label='Estimated Change Point')
    l3, = plt.plot([-1], [-1], color='black', linestyle='dashed', alpha=0.7, label='True Change Point')
    red_patch = matplotlib.patches.Patch(color='red', alpha=0.1, label='95% HPD Interval')

    ax.legend(handles=[l1,l2,red_patch,l3])
    
    
def plot_bg(beta_true, beta_est, beta_mean, beta_mode, beta_median,
            gamma_true, gamma_est, gamma_mean, gamma_mode, gamma_median, T_max=100):
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12,6), layout='constrained')

    for i in range(2):
        for j in range(2):
            ax[i,j].set_xticks(np.arange(5,T_max+1,10))
            ax[i,j].set_xlabel("Day")

    x = np.arange(1, T_max+1)

    ax[0,0].set_title('Smooth Estimator')
    ax[0,0].set_xticks(np.arange(5,T_max+1,10))
    ax[0,0].plot(x, beta_est,   label="Estimated Beta",  color='dodgerblue')
    ax[0,0].plot(x, beta_true,  label="True Beta",       color='dodgerblue', linestyle='dashed', alpha=0.7)
    ax[0,0].plot(x, gamma_est,  label="Estimated Gamma", color='darkorange')
    ax[0,0].plot(x, gamma_true, label="True Gamma",      color='darkorange', linestyle='dashed', alpha=0.7)

    ax[0,1].set_title('Mean')
    ax[0,1].set_xticks(np.arange(5,T_max+1,10))
    ax[0,1].plot(x, beta_mean,  label="Mean Beta",  color='dodgerblue')
    ax[0,1].plot(x, beta_true,  label="True Beta",  color='dodgerblue', linestyle='dashed', alpha=0.7)
    ax[0,1].plot(x, gamma_mean, label="Mean Gamma", color='darkorange')
    ax[0,1].plot(x, gamma_true, label="True Gamma", color='darkorange', linestyle='dashed', alpha=0.7)

    ax[1,0].set_title('Mode')
    ax[1,0].set_xticks(np.arange(5,T_max+1,10))
    ax[1,0].plot(x, beta_mode,  label="Mode Beta",  color='dodgerblue')
    ax[1,0].plot(x, beta_true,  label="True Beta",  color='dodgerblue', linestyle='dashed', alpha=0.7)
    ax[1,0].plot(x, gamma_mode, label="Mode Gamma", color='darkorange')
    ax[1,0].plot(x, gamma_true, label="True Gamma", color='darkorange', linestyle='dashed', alpha=0.7)

    ax[1,1].set_title('Median')
    ax[1,1].set_xticks(np.arange(5,T_max+1,10))
    ax[1,1].plot(x, beta_median,  label="Median Beta",  color='dodgerblue')
    ax[1,1].plot(x, beta_true,    label="True Beta",    color='dodgerblue', linestyle='dashed', alpha=0.7)
    ax[1,1].plot(x, gamma_median, label="Median Gamma", color='darkorange')
    ax[1,1].plot(x, gamma_true,   label="True Gamma",   color='darkorange', linestyle='dashed', alpha=0.7)

    l1, = plt.plot([0], [0], color='dodgerblue')
    l2, = plt.plot([0], [0], color='darkorange')
    l3, = plt.plot([0], [0], color='grey')
    l4, = plt.plot([0], [0], color='grey', linestyle='dashed')

    fig.legend((l1,l2,l3,l4), ('Beta', 'Gamma', 'Estimator', 'True Value'), loc='outside right upper')
    
    
def plot_distrib(chain, var, est, mean, mode, median, bins, times):

    arr = np.array(chain[var])

    nrows = math.ceil(len(times)/2)
    ncols = 2
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10,4*nrows), layout='constrained')

    for i in range(nrows):
        for j in range(ncols):
            n = i*ncols+j
            if n == len(times):
                break
            t = times[n]
            
            ax[i,j].set_title(f'Day {t+1}')
            x, b, _ = ax[i,j].hist(arr[:,t], bins=bins, density=True, color='lightblue')
            x_left = np.min([est[t], b[0]])
            x_right = np.max([est[t], b[-1]])
            x_left = x_left-0.1*(x_right-x_left)
            x_right = x_right+0.1*(x_right-x_left)
            ax[i,j].set_xlim(x_left, x_right)
            height = x[np.argmax(x)] 
            ax[i,j].vlines(x=est[t],    ymin=0, ymax=height, color='darkorange', linestyle='dashed')
            ax[i,j].vlines(x=mean[t],   ymin=0, ymax=height, color='red', linestyle='dashed')
            ax[i,j].vlines(x=mode[t],   ymin=0, ymax=height, color='blue', linestyle='dashed')
            ax[i,j].vlines(x=median[t], ymin=0, ymax=height, color='green', linestyle='dashed')

    l1, = plt.plot([0], [0], color='darkorange', linestyle='dashed')
    l2, = plt.plot([0], [0], color='red', linestyle='dashed')
    l3, = plt.plot([0], [0], color='blue', linestyle='dashed')
    l4, = plt.plot([0], [0], color='green', linestyle='dashed')

    fig.legend((l1,l2,l3,l4), ('Smooth Estimator', 'Mean', 'Mode', 'Median'), loc='outside right upper')