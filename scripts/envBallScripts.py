# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 13:18:48 2025

@author: drgarza
"""

from pathlib import Path
import os

import numpy as np
import pandas as pd


import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

import cobra
from cobra import Model, Reaction
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"

def get_exchange_metabolites(model):
    return set([i.id for i in model.exchanges if i.id != 'EX_biomass(e)'])
    

def apply_env(env_dict, 
             model, 
             upper_bound=1000, 
             make_copy=False):
    
    mc = model.copy() if make_copy else model

    for rxn_id, lb in env_dict.items():
        if mc.reactions.has_id(rxn_id):
            rxn = mc.reactions.get_by_id(rxn_id)
            rxn.lower_bound = -1*abs(lb)
            rxn.upper_bound = upper_bound  
    return mc
    
def apply_env_ball(model, env_ball):
    
    flux_dict = {}

    for env_key, env_sample in env_ball.items():
        m = apply_env(env_sample, model, make_copy=False)
        solution = m.optimize()
        flux_dict[env_key] = solution.fluxes  # pandas Series

    flux_df = pd.DataFrame(flux_dict)
    return flux_df


 
def gen_environment_ball(exchanges,
                         anaerobic=True,
                         fixed_reactions={'EX_h2o(e)': 100},
                         size=1000,
                         total_flux=100,
                         seed=666):
    """
    Generate a dictionary of random environments using a Dirichlet distribution.

    Returns:
        dict of dicts: Keys are environment indices, 
        values are lower-bounds.
    """
    np.random.seed(seed)

    variable_exchanges = [rxn for rxn in exchanges if rxn not in fixed_reactions and rxn != 'EX_o2(e)']
    alpha = np.ones(len(variable_exchanges))
    dirichlet_samples = np.random.dirichlet(alpha, size=size)

    envs = {}

    for i in range(size):
        env_dict = {
            rxn_id: dirichlet_samples[i, j] * total_flux
            for j, rxn_id in enumerate(variable_exchanges)
        }

        # Set oxygen
        env_dict['EX_o2(e)'] = 0.0 if anaerobic else 1000.0

        # Overwrite with fixed reactions
        env_dict.update(fixed_reactions)

        envs[i] = env_dict

    return envs
    

def plot_flux_heatmap(flux_df, 
                      output_path=None, 
                      figsize=(12, 10), 
                      method='average', 
                      metric='correlation'):
   
    # Drop all-zero rows
    data = flux_df.loc[~(flux_df == 0).all(axis=1)]

    # Normalize each row by max absolute value (signed normalization)
    data_norm = data.div(data.abs().max(axis=1), axis=0)

    # Create clustered heatmap
    g = sns.clustermap(data_norm,
                       cmap='vlag',
                       figsize=figsize,
                       method=method,
                       metric=metric,
                       vmin=-1, vmax=1,
                       xticklabels=False, 
                       yticklabels=False,
                       cbar_kws={'label': 'Normalized Flux'},
                       row_cluster=True,
                       col_cluster=True)

    g.ax_heatmap.set_xlabel("Environments")
    g.ax_heatmap.set_ylabel("Reactions")
    g.fig.suptitle("Environment-Driven Flux Clustermap", y=1.02)

    if output_path:
        g.savefig(output_path, dpi=300)

    plt.show()


