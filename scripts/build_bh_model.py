# -*- coding: utf-8 -*-
"""
Created on Fri May 13 13:40:58 2022

@author: u0139894
"""
import cobra
from getModelReactions import *
from cobra import Model

def makeSink(r_id, metabolite):
    sink = Reaction(r_id)
    sink.lower_bound = -1000
    sink.upper_bound = 1000
    sink.add_metabolites({metabolite:-1})
    return sink


modelOutput = 'C:\\Users\\u0139894\\Documents\\GitHub\\metabolic_toy_model\\files\\models\\acetogen_toy_model.xml'

reactions = []
with open('C:\\Users\\u0139894\\Documents\\GitHub\\metabolic_toy_model\\files\\BH_metabolicReactions.txt') as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        if 'rxn' in a[1]:
            reactions.append(a[1])

modelReactions = getModelReactions(reactions)

model = Model("acetogen")
model.add_reactions(modelReactions)

# #add external metabolites
lac_e = creatMetObj('cpd00159', 'e')
co2_e = creatMetObj('cpd00011', 'e')
h2_e = creatMetObj('cpd11640', 'e')
ace_e = creatMetObj('cpd00029', 'e')
for_e = creatMetObj('cpd00047', 'e')
prot_e = creatMetObj('cpd00067', 'e')

#add exchange reactions (using the same fuction for sink)
gluEX = makeSink('EX_cpd00027_e', model.metabolites.cpd00027_e)
lacEX = makeSink('EX_cpd00159_e', lac_e)
co2EX = makeSink('EX_cpd00011_e', co2_e)
h2EX = makeSink('EX_cpd11640_e', h2_e)
aceEX = makeSink('EX_cpd00029_e', ace_e)
forEX = makeSink('EX_cpd00047_e', for_e)
protEX = makeSink('EX_cp00067_e', prot_e)


# #add the RNF reaction

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

# add the R11818 reaction

for_c = model.metabolites.cpd00047_c.copy()
co2_c = model.metabolites.cpd00011_c.copy()
fdh = Reaction('FDH')
fdh.lower_bound = -1000
fdh.upper_bound = 0


fdh.add_metabolites({for_c:-2, nad_c:-1, fdox_c:-2, co2_c:2, nadh_c:1, prot_c:1, fdrd_c:2})

#add H2 transporter
h2T = Reaction('h2T')
h2T.name = 'H2 passive transport'
h2T.lower_bound = -1000
h2T.upper_bound = 1000
h2T.add_metabolites({model.metabolites.cpd11640_c:-1,h2_e:1})

#add CO2 transporter
co2T = Reaction('co2T')
co2T.name = 'co2 passive transport'
co2T.lower_bound = -1000
co2T.upper_bound = 1000

co2T.add_metabolites({model.metabolites.cpd00011_c:-1,co2_e:1})




# #add the objective reaction
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


# #add needed sink reactions
piSink = makeSink('piSink', model.metabolites.cpd00009_c)
h2oSink = makeSink('h2oSink', model.metabolites.cpd00001_c)



#add all new reactions to model

model.add_reactions([gluEX, lacEX, aceEX, co2EX, h2EX, forEX, protEX, rnf, fdh, h2T, co2T, biomass, piSink, h2oSink])



# #############pipe end products to external metabolites#####


#acetate
model.reactions.rxn00225.subtract_metabolites({model.metabolites.cpd00029_c:-1})
model.reactions.rxn00225.add_metabolites({ace_e:-1})
model.reactions.rxn00225.upper_bound=0

#lactate
model.reactions.rxn00499.subtract_metabolites({model.metabolites.cpd00159_c:-1})
model.reactions.rxn00499.add_metabolites({lac_e:-1})
model.reactions.rxn00499.upper_bound = 0

# fix the ferredoxin ids
model.reactions.rxn45849.subtract_metabolites({model.metabolites.cpd27757_c:-2})
model.reactions.rxn45849.subtract_metabolites({model.metabolites.cpd28082_c:2})
model.reactions.rxn45849.add_metabolites({model.metabolites.cpd11621_c:-2})
model.reactions.rxn45849.add_metabolites({model.metabolites.cpd11620_c:2})



####################
model.reactions.rxn05938.lower_bound=-1000
model.reactions.rxn05938.upper_bound=0
####################



#######.

model.reactions.EX_cpd00027_e.lower_bound = -1000
model.reactions.EX_cpd00011_e.lower_bound=-1000
model.reactions.EX_cpd11640_e.lower_bound=-1000

model.reactions.EX_cpd00159_e.lower_bound=0
model.reactions.EX_cpd00029_e.lower_bound=0
model.reactions.EX_cpd00047_e.lower_bound=0
#model.reactions.EX_cp00067_e.lower_bound=0

model.objective = 'biomass'
model.objective.direction = 'max'

cobra.io.write_sbml_model(model, modelOutput)