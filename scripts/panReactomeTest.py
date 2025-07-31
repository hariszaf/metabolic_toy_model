# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 17:55:28 2025

@author: drgarza
"""
import os
import cobra
from pathlib import Path
from build_pan_reactomes import make_panReactome, add_exchange, generate_binary_presence_matrix

root_dir = Path(__file__).resolve().parent.parent
model_folder = root_dir / 'files' / 'models' / 'AGORA' / 'pan_genome'

models = os.listdir(model_folder)

model = cobra.io.read_sbml_model(
    os.path.join(model_folder, 'Bacteroides_thetaiotaomicron_VPI_5482.xml')
)

for reaction in model.reactions:
    if reaction.objective_coefficient == 1.0:
        objective = reaction.copy()

preact = make_panReactome(model_folder, 'Bacteroides_panreactome')

# Path where we save the image and reaction presence absence table
results_folder = root_dir / 'results' / 'pan_reactome'
os.makedirs(results_folder, exist_ok=True)

generate_binary_presence_matrix(
    model_folder,
    preact,
    results_folder
)

add_exchange(preact)
preact.add_reactions([objective])
preact.reactions.get_by_id(objective.id).objective_coefficient = 1.0

solution  = preact.optimize()
print(f"Objective Value: {solution.objective_value: .2f}")

cobra.io.write_sbml_model(preact, os.path.join(results_folder, 'bacteroides_pan_reactome.xml'))
