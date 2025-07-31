# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 13:18:48 2025

@author: drgarza
"""

from pathlib import Path
import os

import cobra
from cobra import Model, Reaction
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"

def get_exchange_metabolites(model):
    return set([i.id for i in model.exchanges if i.id != 'EX_biomass(e)'])
    

def applyEnv(env_dict, 
             model, 
             upper_bound=1000, 
             make_copy=False):
    
    mc = model.copy() if make_copy else model

    for rxn_id, lb in env_dict.items():
        if mc.reactions.has_id(rxn_id):
            rxn = mc.reactions.get_by_id(rxn_id)
            rxn.lower_bound = lb
            rxn.upper_bound = upper_bound  
    return mc
    
def gen_environment_ball(exchanges, 
                         anaerobic = True, 
                         fixed_reactions = {'EX_h2o(e)': 100},
                         size = 1000,
                         seed = 66):
    
 
def gen_environment_ball(exchanges,
                         anaerobic=True,
                         fixed_reactions={'EX_h2o(e)': 100},
                         size=1000,
                         total_flux=100,
                         seed=66):
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
    

root_dir = Path(__file__).resolve().parent.parent
model_folder = root_dir / 'files' / 'models' / 'AGORA' /'no_mucin'


models = [cobra.io.read_sbml_model(model_folder / i) for i in os.listdir(model_folder) if '.xml' in i] 

exchanges = set()

for model in models:
    exchanges = exchanges.union(get_exchange_metabolites(model))