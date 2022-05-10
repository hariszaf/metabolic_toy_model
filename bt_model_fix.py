# -*- coding: utf-8 -*-
"""
Created on Tue May 10 04:59:09 2022

@author: u0139894
"""

import cobra
from cobra import Model, Reaction, Metabolite


def makeSink(r_id, metabolite):
    sink = Reaction(r_id)
    sink.lower_bound = -1000
    sink.upper_bound = 1000
    sink.add_metabolites({metabolite:-1})
    return sink

modelInput = 'bt_toy_model.xml'

modelOutput = 'bt_toy_model_functional.xml'

model = cobra.io.read_sbml_model(modelInput)

##############add external metabolites########

#lactate
lac_e = model.metabolites.cpd00159_c.copy()
lac_e.id = 'cpd00159_e'
lac_e.compartment = 'e'

#formate
for_e = model.metabolites.cpd00047_c.copy()
for_e.id = 'cpd00047_e'
for_e.compartment = 'e'

#acetate
ace_e = model.metabolites.cpd00029_c.copy()
ace_e.id = 'cpd00029_e'
ace_e.compartment = 'e'

#succinate
succ_e = model.metabolites.cpd00036_c.copy()
succ_e.id = 'cpd00036_e'
succ_e.compartment = 'e'


#############################################

#############add exchange reactions############

#lactate
lacExch = Reaction('EX_cpd000159_e')
lacExch.name = 'lactate exchange'
lacExch.lower_bound = -1000
lacExch.upper_bound = 1000

lacExch.add_metabolites({lac_e:-1})


#formate
forExch = Reaction('EX_cpd000047_e')
forExch.name = 'formate exchange'
forExch.lower_bound = -1000
forExch.upper_bound = 1000

forExch.add_metabolites({for_e:-1})

#acetate
aceExch = Reaction('EX_cpd000029_e')
aceExch.name = 'acetate exchange'
aceExch.lower_bound = -1000
aceExch.upper_bound = 1000

aceExch.add_metabolites({ace_e:-1})

#succinate
succExch = Reaction('EX_cpd000036_e')
succExch.name = 'succinate exchange'
succExch.lower_bound = -1000
succExch.upper_bound = 1000

succExch.add_metabolites({succ_e:-1})

#proton
prot_e = model.metabolites.cpd00067_e.copy()
protExch = Reaction('EX_cpd000067_e')
protExch.name = 'proton exchange'
protExch.lower_bound = -1000
protExch.upper_bound = 1000

protExch.add_metabolites({prot_e:-1})



#############add transporters#######
#glucose
prot_c = model.metabolites.cpd00067_c.copy()
glc_c = model.metabolites.cpd00027_c.copy()
glc_e = model.metabolites.cpd00027_e.copy()


glcT = Reaction('rxn05573')
glcT.name = 'D-glucose transport in via proton sympor'
glcT.lower_bound = -1000
glcT.upper_bound = 1000

glcT.add_metabolites({glc_e:-1, prot_e:-1, glc_c:1, prot_c:1})

# #pyruvate
# pyr_c = model.metabolites.cpd00020_c.copy()
# pyr_e = model.metabolites.cpd00020_e.copy()

# pyrT = Reaction('rxn05469')
# pyrT.name = 'pyruvate reversible transport via proton symport'
# pyrT.lower_bound = -1000
# pyrT.upper_bound = 1000

# pyrT.add_metabolites({pyr_e:-1, prot_e:-1, pyr_c:1, prot_c:1})

#add the RNF reaction
nad_c = model.metabolites.cpd00003_c.copy()
nadh_c = model.metabolites.cpd00004_c.copy()
fdox_c = model.metabolites.cpd11621_c.copy()
fdrd_c = model.metabolites.cpd11620_c.copy()
rnf = Reaction('RNF')
rnf.name = 'RNF'
rnf.lower_bound = -1000
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


#add all new reactions to model

model.add_reactions([lacExch, forExch, aceExch, succExch, protExch, glcT, rnf, biomass, piSink, h2oSink])


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


model.objective = 'biomass'

cobra.io.write_sbml_model(model, modelOutput)


