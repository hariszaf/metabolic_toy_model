# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 17:55:28 2025

@author: drgarza
"""

from build_pan_reactomes import *

model_folder = os.path.join(Path(os.getcwd()).parents[0], 
                            'files', 
                            'models', 
                            'AGORA', 
                            'pan_genome')

models = os.listdir(model_folder)

model = cobra.io.read_sbml_model(os.path.join(model_folder, 'Bacteroides_thetaiotaomicron_VPI_5482.xml'))

for reaction in model.reactions:
    if reaction.objective_coefficient==1.0:
        objective = reaction.copy()

preact = make_panReactome(model_folder, 'Bacteroides_panreactome')

#where we save the image and reaction presence absence table
results_folder = os.path.join(Path(os.getcwd()).parents[0], 
                              'files', 
                              'results', 
                              'pan_reactome')
generate_binary_presence_matrix(model_folder, 
                                preact, 
                                results_folder)

add_exchange(preact)
preact.add_reactions([objective])
preact.reactions.get_by_id(objective.id).objective_coefficient = 1.0

solution  = preact.optimize()
print(f"Objective Value: {solution.objective_value: .2f}")

cobra.io.write_sbml_model(preact, os.path.join(results_folder, 'bacteroides_pan_reactome.xml'))
