import os
import pickle
import numpy as np
import pandas as pd

from tqdm import tqdm
from numpy.typing import NDArray
from cobra.core import Reaction, Model
from cobra.flux_analysis import loopless_solution, flux_variability_analysis
from cobra.util.solver import linear_reaction_coefficients
from dingo import MetabolicNetwork, PolytopeSampler


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


def run_pseudo_ll_sampling(
    cobra_model,
    fva_fraction   = 1.0,
    dingo_fraction = 50,
    ess            = 1000,
    sample_limit   = None,
    solver         = None
):
    """
    Performs a pseudo loopless sampling, by setting inner reactions bounds based on loopless FVA 
    """

    # Step 1 — Loopless FVA
    print(f"-- Perform a Flux Variability Analysis with optimal fraction {fva_fraction}.")
    ll_opt_fva = flux_variability_analysis(
        cobra_model, loopless=True, fraction_of_optimum=fva_fraction
    )

    for _, rxn_id in ll_opt_fva.iterrows():
        if rxn_id["minimum"] > rxn_id["maximum"]:
            if abs(rxn_id["minimum"] - rxn_id["maximum"]) < 1e-3:
                rxn_id["minimum"] = rxn_id["maximum"]
            else:
                raise "Loopless FVA returns minimum value rather higher than maximum.."

    # Step 2 — Apply bounds
    print("\n-- Set bounds to inner reactions of the model based on the FVA findings.")
    bounded_model = apply_bounds(cobra_model, ll_opt_fva)

    # Step 3 — Build dingo model
    print(f"\n-- Build a dingo model from the cobra altered one, setting its optimal fraction to {dingo_fraction}.")
    ll_dmodel = MetabolicNetwork.from_cobra_model(bounded_model)
    ll_dmodel.set_opt_percentage(dingo_fraction)
    if solver:
        ll_dmodel.set_solver("gurobi")

    # Step 4 — Sample steady states
    print("\n-- Sample with MMCS")
    sampler = PolytopeSampler(ll_dmodel)
    samples = sampler.generate_steady_states(ess=ess, psrf=True)
    samples = pd.DataFrame(samples, index=[rxn.id for rxn in cobra_model.reactions]).T

    # TODO (Haris Zafeiropoulos, 2025-08-12): make fluxes with super low value equal to 0

    # Step 5 — Loopless optimization on first N samples
    print(
        "\n-- Check whether samples can be considered loopless using the looplsess_solution of cabra."
        "If status is `infeasible` then sample is not ok."
    )
    ll_samples = check_loopy(samples, bounded_model, sample_limit)

    return samples, ll_samples, bounded_model


def apply_bounds(cobra_model, min_max):
    """
    Change lower and upper bound of a cobra model based on a dataframe with reactions as index
    and a `minimum` and `maximum` columns, as the one as cobra.flux_analysis.flux_variability 
    returns.
    """
    cobra_model = cobra_model.copy()

    for rxn, fva_row in zip(cobra_model.reactions, min_max.iterrows()):

        if rxn.id != fva_row[0]:
            print("WE HAVE A MISS-MATCH: ", fva_row.index)
            raise

        if rxn in cobra_model.exchanges:
            print(f"{rxn.id} is among the exchange reactions of the model. Bounds free.")
            continue

        lmin = fva_row[1]["minimum"]
        lmax = fva_row[1]["maximum"]
        if lmin > lmax:
            lmin = lmax

        rxn.lower_bound = lmin
        rxn.upper_bound = lmax

    return cobra_model


def check_loopy(samples, cobra_model, sample_limit = None):
    ll_vs = []
    if sample_limit is None:
        sample_limit = samples.shape[0]
    for _, v in tqdm(samples.iloc[:sample_limit].iterrows(),
                     total=sample_limit,
                     desc=f"Processing first {sample_limit}"):
        try:
            ll_v = loopless_solution(cobra_model, v)  # assumes loopless_solution(model, fluxes)
        except Exception as e:
            print(f"Error in loopless_solution: {e}")
            continue
        if ll_v.status == "optimal":
            ll_vs.append(ll_v)
    return ll_vs
