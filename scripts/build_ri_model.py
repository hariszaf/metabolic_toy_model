# -*- coding: utf-8 -*-
"""
Created on Fri May 13 13:40:58 2022

@author: u0139894
"""

from getModelReactions import *
from cobra import Model

def makeSink(r_id, metabolite):
    sink = Reaction(r_id)
    sink.lower_bound = -1000
    sink.upper_bound = 1000
    sink.add_metabolites({metabolite:-1})
    return sink


modelOutput = 'C:\\Users\\u0139894\\Documents\\GitHub\\metabolic_toy_model\\files\\models\\butyrate_producer_toy_model.xml'
reactions = []
with open('C:\\Users\\u0139894\\Documents\\GitHub\\metabolic_toy_model\\files\\RI_metabolicReactions.txt') as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        if 'rxn' in a[1]:
            reactions.append(a[1])

modelReactions = getModelReactions(reactions)

model = Model("butyrate_producer")
model.add_reactions(modelReactions)

#add external metabolites
lac_e = creatMetObj('cpd00159', 'e')
co2_e = creatMetObj('cpd00011', 'e')
h2_e = creatMetObj('cpd11640', 'e')
ace_e = creatMetObj('cpd00029', 'e')
but_e = creatMetObj('cpd00211', 'e')

#add exchange reactions (using the same fuction for sink)
gluEX = makeSink('EX_cpd00027_e', model.metabolites.cpd00027_e)
lacEX = makeSink('EX_cpd00159_e', lac_e)
co2EX = makeSink('EX_cpd00011_e', co2_e)
h2EX = makeSink('EX_cpd11640_e', h2_e)
aceEX = makeSink('EX_cpd00029_e', ace_e)
butEX = makeSink('EX_cpd00211_e', but_e)
protEX = makeSink('EX_cp00067_e', model.metabolites.cpd00067_e)

#add acetate transporter
aceT = Reaction('aceT')
aceT.name = 'acetate passive transport'
aceT.lower_bound = -1000
aceT.upper_bound = 1000

aceT.add_metabolites({model.metabolites.cpd00029_c:-1,ace_e:1})

#add the RNF reaction
prot_e = model.metabolites.cpd00067_e.copy()
prot_c = model.metabolites.cpd00067_c.copy()
nad_c = model.metabolites.cpd00003_c.copy()
nadh_c = model.metabolites.cpd00004_c.copy()
fdox_c = model.metabolites.cpd11621_c.copy()
fdrd_c = model.metabolites.cpd11620_c.copy()
rnf = Reaction('RNF')
rnf.name = 'RNF'
rnf.lower_bound = -1000
rnf.upper_bound = 1000
rnf.add_metabolites({prot_c:-3,nad_c:-1, fdrd_c:-2, nadh_c:1, prot_e:2, fdox_c:2 })

#add the ferredoxin dependent BTCOADH
crotonoylcoA = model.metabolites.cpd00650_c.copy()
butyrylcoA = model.metabolites.cpd00120_c.copy()
btcoA = Reaction('BTCOAACCOAT')
btcoA.name = 'Butanoyl-CoA:acetate CoA-transferase'
btcoA.lower_bound = 0
btcoA.upper_bound = 1000
btcoA.add_metabolites({crotonoylcoA:-1, fdox_c:-1, prot_c:-2, nadh_c:-2, butyrylcoA:1, fdrd_c:1, nad_c:2})

 




#add the objective reaction
atp_c = model.metabolites.cpd00002_c.copy()
adp_c = model.metabolites.cpd00008_c.copy()
accoA_c = model.metabolites.cpd00022_c.copy()
coA_c = model.metabolites.cpd00010_c.copy()
biomass_c = Metabolite('biomass', compartment='c')
biomass = Reaction('biomass')
biomass.name='Mock biomass function'
biomass.lower_bound=0
biomass.upper_bound=1000
biomass.add_metabolites({atp_c:-3, 
                         accoA_c:-2,
                         nadh_c:-2, 
                         prot_c:-2, 
                         adp_c:3,
                         nad_c:2,
                         coA_c:2
                         })


#add needed sink reactions
piSink = makeSink('piSink', model.metabolites.cpd00009_c)
h2oSink = makeSink('h2oSink', model.metabolites.cpd00001_c)



#add all new reactions to model

model.add_reactions([gluEX, lacEX, aceEX, co2EX, h2EX, butEX, protEX, rnf, btcoA, biomass, piSink, h2oSink])



#############pipe end products to external metabolites#####


#lactate
model.reactions.rxn00499.subtract_metabolites({model.metabolites.cpd00159_c:-1})
model.reactions.rxn00499.add_metabolites({lac_e:-1})
model.reactions.rxn00499.upper_bound = 1000

#co2
model.reactions.rxn05938.subtract_metabolites({model.metabolites.cpd00011_c:-1})
model.reactions.rxn05938.add_metabolites({co2_e:-1})
model.reactions.rxn05938.lower_bound = -1000
model.reactions.rxn05938.upper_bound = 0

#H2 (also fix the ferredoxin ids)
model.reactions.rxn45849.subtract_metabolites({model.metabolites.cpd11640_c:-1})
model.reactions.rxn45849.subtract_metabolites({model.metabolites.cpd27757_c:-2})
model.reactions.rxn45849.subtract_metabolites({model.metabolites.cpd28082_c:2})
model.reactions.rxn45849.add_metabolites({h2_e:-1})
model.reactions.rxn45849.add_metabolites({model.metabolites.cpd11621_c:-2})
model.reactions.rxn45849.add_metabolites({model.metabolites.cpd11620_c:2})
model.reactions.rxn45849.lower_bound = -1000
model.reactions.rxn45849.upper_bound = 0

#unify metabolites with naming problem in the SEEDdb
model.reactions.rxn27735.subtract_metabolites({model.metabolites.cpd02234_c:1})
model.reactions.rxn27735.add_metabolites({model.metabolites.cpd00842_c:1})

#butyrate
model.reactions.rxn00875.subtract_metabolites({model.metabolites.cpd00211_c:1})
model.reactions.rxn00875.add_metabolites({but_e:1})
model.reactions.rxn00875.lower_bound = 0
model.reactions.rxn00875.upper_bound = 1000





#######.

model.reactions.EX_cpd00027_e.lower_bound=-1000
model.reactions.EX_cpd00159_e.lower_bound=0
model.reactions.EX_cpd00011_e.lower_bound=0
model.reactions.EX_cpd11640_e.lower_bound=0
model.reactions.EX_cpd00029_e.lower_bound=0
model.reactions.EX_cpd00211_e.lower_bound=0

model.reactions.rxn02167.lower_bound=0

model.objective = 'biomass'
model.objective.direction = 'max'

cobra.io.write_sbml_model(model, modelOutput)