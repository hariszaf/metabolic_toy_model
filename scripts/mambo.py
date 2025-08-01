from pathlib import Path
import os

import numpy as np
import pandas as pd


import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np
import scipy.stats as sts

import cobra
from cobra import Model, Reaction
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"

def get_exchange_metabolites(model):
    return set([i.id for i in model.exchanges if i.id != 'EX_biomass(e)'])

def apply_environment(mdl, media_dict):
    for i in media_dict:
        if mdl.reactions.has_id(i):
            mdl.reactions.get_by_id(i).lower_bound=-media_dict[i]

    sol_m=mdl.optimize()
    return sol_m.objective_value

def getRMSE(x,y):

    return np.sqrt( np.mean( (x-y)**2) )

def MetropolisSA(new_value, current_value, rabs, T):

    if None in new_value:
        return (False, 0, 0, 1)

    candidate_prob = np.exp(-getRMSE(current_value, new_value) / T)


    if np.random.uniform() < candidate_prob:
        return (True, candidate_prob) #True means the statistics will be accepted.
    else:
        return (False, candidate_prob)#Statistics rejected.


def Metropolis(new_value, current_value, rabs):

    if None in new_value:
        return (False, 0, 0, 1)
    optimal = np.arcsinh(0.999)
    standard_error = 0.45#(1/(len(new_value)**0.5))
    candidate_prob = sts.norm.pdf(np.arctanh(sts.pearsonr(new_value, rabs)[0]), loc=optimal,scale=standard_error)

    current_prob = sts.norm.pdf(np.arctanh(sts.pearsonr(current_value, rabs)[0]), loc=optimal,scale=standard_error)

    statistic = candidate_prob/current_prob

    if np.random.uniform() < statistic:
        return (True, statistic, candidate_prob, current_prob) #True means the statistics will be accepted.
    else:
        return (False, statistic, candidate_prob, current_prob)#Statistics rejected.


def current_solution(modelList, media):
    sol = np.array([apply_environment(i, media) for i in modelList])
    return sol


def MCMC(media, modelList, rab, delta = 0.5):
    m2 = media.copy()
    ch = np.random.choice(list(media))
    m2[ch]= max(0, m2[ch] + np.random.uniform(low=-delta, high=delta))
    sol_current = current_solution(modelList, media)
    sol_candidate = current_solution(modelList, m2)

    met = Metropolis(sol_candidate, sol_current, rab)


    print(met)

    if met[0]:
        return (sol_candidate, m2)

    else:
        return (sol_current, media)

def bunching(vec):
    p0 = vec[0]

    for i in range(1, len(vec)):
        p0 = (p0 + vec[i])/2

    return p0


def plot_mambo_results(result, cor_sorter, output_path=None, figsize=(12, 6), cmap='coolwarm'):
    
    fig, ax1 = plt.subplots(figsize=figsize)

    # Plot pcolormesh of media matrix
    c = ax1.pcolormesh(result, cmap=cmap, shading='auto')
    ax1.set_ylabel("Metabolites (media components)")
    ax1.set_xlabel("MCMC Samples (sorted by correlation)")
    ax1.set_title("Sorted MAMBO Media vs. Correlation with Composition")

    # Add colorbar for media intensities
    cbar = plt.colorbar(c, ax=ax1)
    cbar.set_label("Media concentration")

    # Plot overlayed correlation line
    ax2 = ax1.twinx()
    ax2.plot(cor_sorter, color='black', linewidth=2.5, label='Pearson correlation with target')
    ax2.set_ylabel("")
    ax2.set_ylim(-1, 1)

    # Optional: legend for the line
    ax2.legend(loc='upper right')

    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300)
    plt.show()
