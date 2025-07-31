# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 15:08:01 2025

@author: drgarza
"""

from envBallScripts import *

root_dir = Path(__file__).resolve().parent.parent
model_folder = root_dir / 'files' / 'models' / 'AGORA' /'no_mucin'


model = cobra.io.read_sbml_model(model_folder / 'Bacteroides_thetaiotaomicron_VPI_5482.xml')

exchanges = get_exchange_metabolites(model)
    
envBall = gen_environment_ball(exchanges,
                         anaerobic=True,
                         fixed_reactions={'EX_h2o(e)': 100},
                         size=1000,
                         total_flux=100,
                         seed=666)

solutions = apply_env_ball(model, envBall)

results_folder = root_dir / 'results' / 'env_ball'
os.makedirs(results_folder, exist_ok=True)

outputPath = results_folder / 'env_ball_reactions_cluster_bt.png'

plot_flux_heatmap(solutions, outputPath)