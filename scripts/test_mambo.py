# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 13:36:18 2025

@author: drgarza
"""

from mambo import *
from envBallScripts import *



root_dir = Path(__file__).resolve().parent.parent
model_path = root_dir / 'files' / 'models' / 'agora' / 'no_mucin'


acetogen = cobra.io.read_sbml_model(model_path / 'Blautia_hydrogenotrophica_DSM_10507.xml')

sugar_fermenter = cobra.io.read_sbml_model(model_path / 'Bacteroides_thetaiotaomicron_VPI_5482.xml')

butyrate_producer = cobra.io.read_sbml_model(model_path / 'Roseburia_intestinalis_L1_82.xml')


exchanges = get_exchange_metabolites(acetogen)

exchanges = exchanges.union(get_exchange_metabolites(sugar_fermenter))

exchanges = exchanges.union(get_exchange_metabolites(butyrate_producer))


media = {i: 1 for  i in exchanges}

apply_environment(acetogen, media)
apply_environment(sugar_fermenter, media)
apply_environment(butyrate_producer, media)


acetogen.optimize()
sugar_fermenter.optimize()
butyrate_producer.optimize()


modelList = [sugar_fermenter, butyrate_producer, acetogen]

# #####Composition Vector####
composition = np.array([10, 1, 0.1])

medias = [np.array(list(media.values()))]
solutions = [current_solution(modelList, media)]

for i in range(10000):#should be much larger

    print(i)
    solution, media = MCMC(media, modelList, composition)

    if (i>10):#should be much larger

        medias.append(np.array(list(media.values())))
        solutions.append(solution)
        
        


medias= np.array(medias)
mediasM = medias.copy()

medias= medias.T

maxMedias = np.max(mediasM, axis=1)
mediasM = np.array([mediasM[i]/maxMedias[i] for i in range(len(maxMedias))]).T

mSols = solutions.copy()
solutions = np.array(solutions).T
maxSolutions = np.max(mSols, axis=1)
mSols = np.array([mSols[i]/maxSolutions[i] for i in range(len(maxSolutions))]).T


cor = np.array([sts.pearsonr(i, composition)[0] for i in solutions.T])

print(cor)


avM = np.median(medias.T[cor>0.6], axis=0)
avM = (avM/max(avM))*10
m = {list(media.keys())[i]: avM[i] for i in range(len(avM))}

print(f"composition was {composition} \t MAMBO solution was: {current_solution(modelList, m)}")