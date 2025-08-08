# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 18:23:12 2025

@author: drgarza
"""

import cobra
import numpy as np
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt

from .envBallScripts import apply_env


def remove_reaction(orig_obj, model, reaction_id):
    """
    Remove a reaction from the model and check if the objective value increases.
    If it does, return 0, otherwise return 1.
    """
    up_orig = model.reactions.get_by_id(reaction_id).upper_bound
    lo_orig = model.reactions.get_by_id(reaction_id).lower_bound

    model.reactions.get_by_id(reaction_id).upper_bound = 0
    model.reactions.get_by_id(reaction_id).lower_bound = 0

    new_obj = np.round(model.slim_optimize(),decimals=5)

    if new_obj > orig_obj:
        return 0.0

    else:
        model.reactions.get_by_id(reaction_id).upper_bound = up_orig
        model.reactions.get_by_id(reaction_id).lower_bound = lo_orig
        return 1.0


def get_panEFM(map_d, ordered_reactions, modelPath, cuttof):
    '''
    Parameters
    ----------
    map_d : dict
        assures a consistent order of reactions
    ordered_reactions : list
        an arbitrary reaction order to be removed one at a time
    model : cobra.model
    cuttof : float
        values where to consider the absence of growth

    Returns
    -------
    essential_r : numpy.array
        a binary vector ordered according to the map_d of irreducible reactions
    '''

    mod = cobra.io.read_sbml_model(modelPath)

    essential_r = np.zeros(len(map_d))
    orig_obj    = np.round(cuttof * mod.slim_optimize(), decimals=5)

    for i in ordered_reactions:
        c                     = remove_reaction(orig_obj, mod,i)
        essential_r[map_d[i]] = c

    return essential_r


def get_panEFM_dist(modelPath, reactions, env_dict, max_it=1000):

    model = cobra.io.read_sbml_model(modelPath)

    rc    = reactions[:]
    map_d = {reactions[i]: i for i in range(len(reactions))}

    apply_env(
        env_dict,
        model,
        upper_bound = 1000,
        make_copy   = False
    )

    all_iters = []

    for i in tqdm(range(max_it)):
        np.random.shuffle(rc)
        panEFM = get_panEFM(map_d, rc, modelPath, 0.01)
        all_iters.append(panEFM)
    return all_iters


def plot_reaction_freq_heatmap(frequency_df, output_path=None, figsize=(12, 8), cmap="viridis"):

    # sort reactions by mean frequency
    sorted_cols = frequency_df.mean(axis=0).sort_values(ascending=False).index
    data = frequency_df[sorted_cols]

    plt.figure(figsize=figsize)
    ax = sns.heatmap(
        data,
        cmap=cmap,
        cbar_kws={'label': 'Reaction frequency'},
        xticklabels=False,
        yticklabels=True
    )

    ax.set_xlabel("Reactions (sorted by frequency)")
    ax.set_ylabel("Environments")
    ax.set_title("Pan-EFM Reaction Frequencies")
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300)
    plt.show()
