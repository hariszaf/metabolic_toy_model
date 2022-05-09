#!/usr/bin/env python3
from cobra import Model, Reaction, Metabolite
from openpyxl import load_workbook
import json, sys, re

# Load reactions from the ModelSEED database (reference)
f = open("ModelSEED_db_reactions.json", "r")
modelSEED_reactions = json.load(f)

# Load compounds from the ModelSEED database (reference)
g = open("ModelSEED_db_compounds.json", "r")
modelSEED_compounds = json.load(g)

# Function to clean up a list 
def remove_blanks_and_pluses(a_list):
    for index, entry in enumerate(a_list):
        repl = entry
        if entry.count("+") > 1:
            repl = entry.split("+")[0]
        repl = repl.strip()
        if "] +" in repl:
            repl = repl[:-2]
        a_list[index] = repl
    return a_list

# --------------------
# Part A: ModelSEED ids - names dictionary 
# --------------------

# Parse ModelSEED database reactions to map compounds with their corresponding ids 
compounds_to_ids = {}
for reaction in modelSEED_reactions:

    definition = reaction["definition"]
    equation   = reaction["equation"]

    definition = re.split('<=>|<=|=>', definition)

    def_comp = remove_blanks_and_pluses(list(filter(None, re.compile("\(\d\)").split(definition[0]))))
    def_prod = remove_blanks_and_pluses(list(filter(None, re.compile("\(\d\)").split(definition[1]))))

    compounds_names = []
    for x in list(filter(None, def_comp + def_prod)):
        compounds_names.append(x.split("[")[0])

    compounds_ids   = []
    for entry in equation.split(" "):
        if "cpd" in entry:
            compounds_ids.append(entry.split("[")[0])

    for k, v in zip(compounds_names, compounds_ids):
        compounds_to_ids[k] = v

with open('comp_names_to_modelSEED_ids.json', 'w') as outfile:
    json.dump(compounds_to_ids, outfile)



# --------------------
# Part B
# --------------------


# Init model 
model = Model('toy_bt')

# Set cases of extracellular 


# Read toy model's reactions 
my_reactions_excel = load_workbook(filename = 'metabolicReactions.xlsx')
bt_reactions_sheet = my_reactions_excel["bt"]

# Reading sheet per row
bt_reactions_with_rxn_id = []
manual_reactions         = []
set_of_metabolites       = set()

for row in bt_reactions_sheet.iter_rows(min_row=1, max_col=5, values_only=True):

    reaction_name = row[2]
    if row[1] != None:
        for entry in modelSEED_reactions:
            """
            example of an entry:
            {'abbreviation': 'GLCt2', 'abstract_reaction': None, 
            'aliases': ['BiGG: GLCt2; GLCt2pp', 'MetaCyc: RXN0-7077', 'iAF1260: GLCt2pp', 'iAG612: rbs_312; rbs_313', 'iAO358: rll_538', 'Name: D-glucose transport in via proton symport; D-glucose:proton symport; D-glucosetransportinviaprotonsymport'], 
            'code': '(1) cpd00027[1] <=> (1) cpd00027[0]', 
            'compound_ids': 'cpd00027;cpd00067', 
            'definition': '(1) D-Glucose[1] + (1) H+[1] <=> (1) D-Glucose[0] + (1) H+[0]', 
            'deltag': 0.0, 'deltagerr': 2.18, 'direction': '=', 'ec_numbers': None, 
            'equation': '(1) cpd00027[1] + (1) cpd00067[1] <=> (1) cpd00027[0] + (1) cpd00067[0]', 
            'id': 'rxn05573', 'is_obsolete': 0, 'is_transport': 1, 'linked_reaction': 'rxn08616', 'name': 'D-glucose transport in via proton symport', 
            'notes': ['GCC', 'EQC'], 'pathways': ['MetaCyc: Alcohol-Biosynthesis (Fermentation to Alcohols); Energy-Metabolism (Generation of Precursor Metabolite and Energy); Fermentation (); PWY-7385 (1,3-propanediol biosynthesis (engineered))'], 
            'reversibility': '=', 'source': 'Primary Database', 'status': 'OK', 
            'stoichiometry': '-1:cpd00027:1:0:"D-Glucose";-1:cpd00067:1:0:"H+";1:cpd00027:0:0:"D-Glucose";1:cpd00067:0:0:"H+"'}
            """
            if entry["name"] == reaction_name:
                entry["upper_bound"] = 1000
                entry["lower_bound"] = -1000

                for participant in entry["stoichiometry"].split(";"):
                    metabolite    = participant.split(":")[1]
                    if metabolite in set_of_metabolites: 
                        continue
                    for term in modelSEED_compounds:
                        if term["id"] == metabolite:
                            set_of_metabolites.add(metabolite)
                            if row[4] == metabolite:
                                comp = "e"
                            else:
                                comp = "c"
                            model_metabolite = Metabolite(term["id"], name = term["name"], formula = term["formula"], compartment = comp)
                            model.add_metabolites([model_metabolite])

                bt_reactions_with_rxn_id.append(entry)
    else:
        manual_reactions.append(row[3])

# Parse manual reactions
manual_compounds = []

for reaction in manual_reactions:
    parts = re.split('<=>|<=|=>', reaction)
    compounds = []
    for c in parts:
        tmp = c.split("+")
        for t in tmp: 
            compounds.append(t)
    compounds = list(filter(None, compounds))
    tmp = []
    for i in compounds: 
        z = re.compile("\(\d\)").split(i)
        for j in z: 
            tmp.append(j.strip())

    clean_compounds = list(filter(None, tmp))
    for i in clean_compounds:
        if i not in manual_compounds:
            manual_compounds.append(i)

# Init model using the reactions that correspond to ModelSEED ids
for case in bt_reactions_with_rxn_id:

    reaction = Reaction(case["id"])
    model.add_reactions([reaction])
    reaction.name = case["name"]
    reaction.lower_bound = case["lower_bound"]
    reaction.upper_bound = case["upper_bound"]

    for participant in case["stoichiometry"].split(";"):
        metabolite    = participant.split(":")[1]
        stoichiometry = participant.split(":")[0]
        reaction.add_metabolites({metabolite : int(stoichiometry)})


# Parse manual reactions to link their compounds to the ModelSEED ones
for compound in manual_compounds: 
    for modelSEED_compound in modelSEED_compounds: 
        if modelSEED_compound["name"] == compound: 
            model_metabolite = Metabolite(term["id"], name = term["name"], formula = term["formula"])
            model.add_metabolites([model_metabolite])
            break

# Add the biomass and other non ModelSEED reactions
for reaction in manual_reactions:
    print(reaction)
sys.exit(0)
# model.objective = 'biomass'

# ---------------------
# Model validation
# ---------------------

# Get ModelSEED ids for those compounds on the manual reactions - if available
for compound in manual_compounds:
    try:
        print(compounds_to_ids[compound])
    except:
        print(compound, " not in ModelSEED")

import tempfile
from pprint import pprint
from cobra.io import write_sbml_model, validate_sbml_model
with tempfile.NamedTemporaryFile(suffix = '.xml') as f_sbml:
    write_sbml_model(model, filename = f_sbml.name)
    report = validate_sbml_model(filename = f_sbml.name)

pprint(report)





