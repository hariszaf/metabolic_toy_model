import os
import pickle
import numpy as np
import pandas as pd

from tqdm import tqdm
from numpy.typing import NDArray
from cobra.core import Reaction, Model
from cobra.flux_analysis import loopless_solution
from cobra.util.solver import linear_reaction_coefficients


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


# Save samples to pickle file
def dump_samples(samples, filename):
    with open(filename, "wb") as f:
        pickle.dump(samples, f)


# Load pickle file
def load_samples(filename):
    with open(filename, "rb") as f:
        df = pickle.load(f)
    return df


def get_loopless_solutions_from_samples(
    samples: NDArray[np.float64],
    cobra_model: Model
) -> NDArray[np.float64]:
    """
    Function that calculates the `loopless_solution` from each sample and saves results into a new numpy 2D array

    Keyword arguments:
    samples (NDArray[np.float64]) -- Numpy 2D array of the samples
    cobra_model (model) -- cobra model object

    Returns:
    samples_loopless_solutions (NDArray[np.float64]) -- Numpy 2D array of the samples after application of `loopless_solution`
    """

    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    samples = samples.copy()

    if samples.shape[0] == samples.shape[1]:
        raise ValueError("Samples array provided has equal rows and columns dimensions. Please change the number of samples")

    # if provided sampling dataset has reactions as rows ==> transpose
    if samples.shape[0] == len(cobra_reactions_str):
        samples = samples.T

    loopless_solutions_pandas_series_default = []
    for i in tqdm(range(samples.shape[0]), desc="Processing samples"):
        sample = samples[i]
        sample_reactions_dictionary = {k:v for k,v in zip(cobra_reactions_str, sample)}
        try:
            loopless_sample = loopless_solution(model=cobra_model, fluxes=sample_reactions_dictionary)
        except Exception:
            pass
        loopless_solutions_pandas_series_default.append(loopless_sample.fluxes)

    df = pd.concat(loopless_solutions_pandas_series_default, axis=1)
    samples_loopless_solutions = df.to_numpy()

    return samples_loopless_solutions


def get_objective_functions(cobra_model: Model):

    objectives_dict         = linear_reaction_coefficients(cobra_model)
    objective_functions     = list(objectives_dict.keys())
    objective_functions_ids = [rxn.id for rxn in objective_functions]

    return objective_functions_ids
