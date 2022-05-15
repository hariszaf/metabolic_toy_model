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


modelOutput = 'C:\\Users\\u0139894\\Documents\\GitHub\\metabolic_toy_model\\files\\models\\sugar_fermenter_toy_model.xml'
reactions = []
with open('C:\\Users\\u0139894\\Documents\\GitHub\\metabolic_toy_model\\files\\BT_metabolicReactions.txt') as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        if 'rxn' in a[1]:
            reactions.append(a[1])

modelReactions = getModelReactions(reactions)

model = Model("sugar_fermenter")
model.add_reactions(modelReactions)

##############add external metabolites########

#lactate
lac_e = creatMetObj('cpd00159', 'e')

#formate
for_e = creatMetObj('cpd00047', 'e')

#acetate
ace_e = creatMetObj('cpd00029', 'e')

#succinate
succ_e = creatMetObj('cpd00036', 'e')

#############################################


#############add exchange reactions############
#glucose
gluEX = makeSink('EX_cpd00027_e', model.metabolites.cpd00027_e)

#lactate
lacEX = makeSink('EX_cpd00159_e', lac_e)

#formate
forEX = makeSink('EX_cpd00047_e', for_e)

#acetate
aceEX = makeSink('EX_cpd00029_e', ace_e)

#succinate
succEX = makeSink('EX_cpd00036_e', succ_e)

#proton
prot_e = model.metabolites.cpd00067_e.copy()
protEX = makeSink('EX_cpd00067_e', prot_e)

#add the RNF reaction
prot_e = model.metabolites.cpd00067_e.copy()
prot_c = model.metabolites.cpd00067_c.copy()
nad_c = model.metabolites.cpd00003_c.copy()
nadh_c = model.metabolites.cpd00004_c.copy()
fdox_c = model.metabolites.cpd11621_c.copy()
fdrd_c = model.metabolites.cpd11620_c.copy()
rnf = Reaction('RNF')
rnf.name = 'RNF'
rnf.lower_bound = 0
rnf.upper_bound = 1000
rnf.add_metabolites({prot_c:-3,nad_c:-1, fdrd_c:-2, nadh_c:1, prot_e:2, fdox_c:2 })


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
protSink = makeSink('protSink', model.metabolites.cpd00067_c)
atpSink = makeSink('atpSink', model.metabolites.cpd00002_c)
adpSink = makeSink('adpSink', model.metabolites.cpd00008_c)
nadSink = makeSink('nadSink', model.metabolites.cpd00003_c)
nadhSink = makeSink('nadhSink', model.metabolites.cpd00004_c)
xSink = makeSink('xSink', model.metabolites.cpd00032_c)



#add all new reactions to model

model.add_reactions([gluEX, lacEX, forEX, aceEX, succEX, protEX, rnf, biomass, piSink, h2oSink, xSink, nadSink, nadhSink])



#############pipe end products to external metabolites#####
#acetate
model.reactions.rxn00225.subtract_metabolites({model.metabolites.cpd00029_c:-1})
model.reactions.rxn00225.add_metabolites({ace_e:-1})
model.reactions.rxn00225.upper_bound=0

#lactate
model.reactions.rxn00499.subtract_metabolites({model.metabolites.cpd00159_c:-1})
model.reactions.rxn00499.add_metabolites({lac_e:-1})
model.reactions.rxn00499.upper_bound = 0


#formate
model.reactions.rxn00157.subtract_metabolites({model.metabolites.cpd00047_c:-1})
model.reactions.rxn00157.add_metabolites({for_e:-1})
model.reactions.rxn00157.upper_bound = 0

#succinate
model.reactions.rxn00284.subtract_metabolites({model.metabolites.cpd00036_c:-1})
model.reactions.rxn00284.add_metabolites({succ_e:-1})
model.reactions.rxn00284.upper_bound = 0
#####################################################################
model.reactions.rxn13974.lower_bound=-1000

####################################################################

model.reactions.EX_cpd00027_e.lower_bound=-1000
model.reactions.EX_cpd00159_e.lower_bound=0
model.reactions.EX_cpd00047_e.lower_bound=0
model.reactions.EX_cpd00029_e.lower_bound=0
model.reactions.EX_cpd00036_e.lower_bound=0
#model.reactions.EX_cpd00067_e.lower_bound = 0



model.objective = 'biomass'
model.objective.direction = 'max'

cobra.io.write_sbml_model(model, modelOutput)