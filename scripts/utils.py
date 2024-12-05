import os
from cobra.core import Reaction

def get_root_dir_from_script():
    # Get the absolute path of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Walk up the directory tree until we find 'metabolic_toy_model'
    while script_dir:
        if 'metabolic_toy_model' in os.path.basename(script_dir):
            return script_dir
        script_dir = os.path.dirname(script_dir)
    return None

def makeSink(r_id, metabolite):
    sink = Reaction(r_id)
    sink.lower_bound = -1000
    sink.upper_bound = 1000
    sink.add_metabolites({metabolite:-1})
    return sink


def get_reactions_from_tsv(path_to_file):
    reactions = []
    with open(path_to_file) as f:
        f.readline()
        for line in f:
            a = line.strip().split('\t')
            if 'rxn' in a[1]:
                reactions.append(a[1])
    return reactions


def apply_environment(model, new_media):

    model.medium = new_media
    return model
