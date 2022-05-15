# -*- coding: utf-8 -*-
"""
Created on Fri May 13 09:54:48 2022

@author: u0139894
"""

import cobra
from cobra import Model, Reaction, Metabolite

from MSEED_compounds import Compounds
from MSEED_reactions import Reactions

#Use script from ModelSEED biochemistry to parse all metabolite/reaction info
compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
compounds_helper.saveCompounds(compounds_dict)
compounds_aliases_dict = compounds_helper.loadMSAliases()
compounds_helper.saveAliases(compounds_aliases_dict)
# compounds_names_dict = compounds_helper.loadNames()
# compounds_helper.saveNames(compounds_names_dict)


reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
reactions_aliases_dict = reactions_helper.loadMSAliases()
reactions_helper.saveAliases(reactions_aliases_dict)


def creatMetObj(seedID, compartment):
    m = compounds_dict[seedID]
    if compartment=='e':
        cid = seedID + '_e'
    else:
        cid = seedID + '_c'
    
    met = Metabolite(cid)
    met.formula = compounds_dict[seedID]['formula']
    met.name = compounds_dict[seedID]['name']
    met.compartment = compartment
    return met

def getMetabolites(reactionList):
    modelMets = {}
    for reac in reactionList:
        s = reactions_dict[reac]['stoichiometry'].split(';')
        for met in s:
            m = met.split(':')
            if m[2] =='1':
                cpd = creatMetObj(m[1], 'e')
                modelMets[cpd.id] = cpd.copy()
            else:
                cpd = creatMetObj(m[1], 'c')
                modelMets[cpd.id] = cpd.copy()
    return modelMets         
    
def getModelReactions(reactionsList):
    reactions = []
    modelMets = getMetabolites(reactionsList)
    for reac in reactionsList:
        R = Reaction(reac)
        R.name =reactions_dict[reac]['name']
        s = reactions_dict[reac]['stoichiometry'].split(';')
        for met in s:
            m = met.split(':')
            if m[2] == '1':
                R.add_metabolites({modelMets[m[1]+'_e']:float(m[0])})
            else:
                R.add_metabolites({modelMets[m[1]+'_c']:float(m[0])})
        if reactions_dict[reac]['reversibility'] == '=':
            R.lower_bound = -1000
        R.upper_bound = 1000
        
        reactions.append(R.copy())
    
    return reactions
        
        