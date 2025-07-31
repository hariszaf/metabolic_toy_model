# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 18:43:03 2025

@author: drgarza
"""

from get_panEFMs import *

root_dir = Path(__file__).resolve().parent.parent
model_path = root_dir / 'results' / 'pan_reactome' / 'bacteroides_pan_reactome.xml'


model = cobra.io.read_sbml_model(model_path)

exchanges = get_exchange_metabolites(model)
    
envBall = gen_environment_ball(exchanges,
                         anaerobic=True,
                         fixed_reactions={'EX_h2o(e)': 100},
                         size=250,
                         total_flux=100,
                         seed=666)




reactions = [i.id for i in model.reactions if i.objective_coefficient == 0]



pan_efms = {}

for i in envBall:
    print(f"Environment Simulation: {i}")
    
    pan_efms[i] = get_panEFM_dist(model_path, reactions, envBall[i], max_it=100)