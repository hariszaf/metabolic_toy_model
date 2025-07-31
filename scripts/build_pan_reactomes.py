# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 14:51:50 2025

@author: drgarza
"""

from pathlib import Path
import os

import cobra
from cobra import Model, Reaction
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list

def generate_binary_presence_matrix(model_folder, panreactome_model, output_path=None):
    """
    Generates a binary matrix comparing which models contain each reaction from 
    the panreactome model.
    Rows (models) are clustered by similarity, columns are sorted by frequency.
    """
    model_files = [f for f in os.listdir(model_folder) if f.endswith('.xml')]
    panreactome_rxns = [r.id for r in panreactome_model.reactions]

    # Initialize matrix
    binary_matrix = pd.DataFrame(0, index=model_files, columns=panreactome_rxns)

    for file in model_files:
        path = os.path.join(model_folder, file)
        model = cobra.io.read_sbml_model(path)  # <-- kept exactly as in your code
        model_rxns = {rxn.id for rxn in model.reactions}
        shared_rxns = model_rxns.intersection(panreactome_rxns)
        binary_matrix.loc[file, list(shared_rxns)] = 1

    # Sort columns by frequency
    col_order = binary_matrix.sum().sort_values(ascending=False).index
    binary_matrix = binary_matrix[col_order]

    # Cluster rows by similarity (Hamming distance)
    if len(binary_matrix) > 1:  # only if more than one model
        distance_matrix = pdist(binary_matrix.values, metric="hamming")
        linkage_matrix = linkage(distance_matrix, method="average")
        row_order = leaves_list(linkage_matrix)
        binary_matrix = binary_matrix.iloc[row_order]

    # Plot
    plt.figure(figsize=(12, 0.25 * len(binary_matrix)))
    sns.heatmap(binary_matrix, cmap="binary", cbar=False)
    plt.xlabel("Reactions (sorted by frequency)")
    plt.ylabel("Models")
    plt.title("Reaction Presence/Absence (clustered)")
    plt.tight_layout()

    if output_path:
        binary_matrix.to_csv(output_path + "_binary_matrix.csv")
        plt.savefig(output_path + "_heatmap.png", dpi=300)

    plt.show()
    return binary_matrix



def clone_model_without_transp_and_obj(model):
    mc = model.copy()
    
    reactions_to_remove = [reac for reac in mc.exchanges]
    for reaction in mc.reactions:
        if reaction.objective_coefficient == 1.0:
            reactions_to_remove.append(reaction)
    mc.remove_reactions(reactions_to_remove, remove_orphans=1)
    
    mc.optimize()
    
    return mc


def make_panReactome(model_folder, modelID):
    
    files = [i for i in os.listdir(model_folder) if 'xml' in i]
    
    #open the first model
    model = cobra.io.read_sbml_model(os.path.join(model_folder, files[0]))
    #join all non exchange reactopms from models
    mc= clone_model_without_transp_and_obj(model)
    mc.id = modelID
    mc.id 
    for i in files:
        mod= cobra.io.read_sbml_model(os.path.join(model_folder, i))
        clone_mod = clone_model_without_transp_and_obj(mod)
        for reaction in clone_mod.reactions:
            if not mc.reactions.has_id(reaction.id):
                mc.add_reactions([reaction.copy()])
                mc.repair() 
    return mc

def add_exchange(model):
    exchange_mets = [metab for metab in model.metabolites if metab.compartment == 'e']
    for metab in exchange_mets: 
        react = Reaction('EX_' + metab.id.replace('[e]','') + '(e)')
        react.name = 'export of ' + metab.name
        react.lower_bound = -1000.  # This is the default
        react.upper_bound = 1000.  # This is the default
        react.add_metabolites({metab: -1.0})
        model.add_reactions([react])
                
        model.repair()
        model.optimize()

