import cobra
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats
from cobra import Model
from typing import Tuple, List
from numpy.typing import NDArray
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests


# Apply PCA analysis on a samples dataset
def pca_samples(samples):
    """
    samples (pd.DataFrame) -- a df with samples as rows and reactions as column names
    """

    # Perform PCA
    pca        = PCA(n_components=10)  # You can change the number of components if needed
    pca_result = pca.fit_transform(samples)

    # Explained variance ratio
    explained_variance = pca.explained_variance_ratio_
    print(f"Explained Variance Ratio: {explained_variance}")

    # Plot the first two principal components
    plt.figure(figsize=(8, 6))
    plt.scatter(pca_result[:, 0], pca_result[:, 1])
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('PCA: First Two Principal Components')
    plt.grid(True)
    plt.show()

    # Optional: Cumulative explained variance (for how much variance is explained by the first n components)
    cumulative_variance = np.cumsum(explained_variance)
    plt.figure(figsize=(8, 6))
    plt.plot(range(1, len(cumulative_variance) + 1), cumulative_variance, marker='o')
    plt.xlabel('Number of Components')
    plt.ylabel('Cumulative Explained Variance')
    plt.title('Cumulative Explained Variance by PCA Components')
    plt.grid(True)
    plt.show()


# Function from dingo-stats:
# github.com/SotirisTouliopoulos/dingo-stats/blob/master/src/distributions_comparison_utils.py
def significantly_altered_reactions(
    conditions          : List[NDArray[np.float64]] = [],
    selected_comparisons: Tuple                     = [(0, 1)],
    cobra_model         : Model                     = None,
    p_value_cutoff      : float                     = 0.05,
    fold_change_cutoff  : float                     = 1.2,
    std_cutoff          : float                     = 1e-2
) -> Tuple[List[str], List[str]]:
    """
    Function that takes as input at least 2 flux sampling conditions to compare and identifies significantly altered reactions
    It performs a Kolmogorov-Smirnov (KS) non-parametric test and corrects p-value for multiple comparisons.
    It additionally calculates a fold change that together with the p-value classifies reactions as significantly altered or not. 

    Keyword arguments:
    conditions (List) -- List of different flux sampling arrays
    selected_comparisons (List) -- List showing which conditions to compare (useful when comparing more than 2 sampling arrays)
    cobra_model (Model) -- cobra model object
    p_value_cutoff (float) -- cutoff for p value from KS test to consider 2 distributions significantly different
    fold_change_cutoff (float) -- cutoff for fold-change to consider 2 distributions significantly different
    std_cutoff (float) -- cutoff to ensure distributions that are compared are not fixed to a certain value

    Returns:
    Tuple[List, List]
        significant_diff_reactions (List) -- List containing reactions significantly altered between the conditions
        not_significant_diff_reactions (List) -- List containing reactions not significantly altered between the conditions
        pval_dict (Dict) -- Dictionary mapping reaction ID to corrected p-value.
        fold_change_dict (Dict) -- Dictionary mapping reaction ID to fold-change value.
    """

    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]

    # if provided sampling dataset has reactions as cols ==> transpose
    for condition in conditions:
        if condition.shape[0] == condition.shape[1]:
            raise ValueError("Samples array provided has equal rows and columns dimensions. Please change the number of samples")

        if condition.shape[1] == len(cobra_reactions_str):
            condition = condition.T

    # Store p-values for each row for KS test
    p_values = {row: [] for row in range(len(cobra_reactions_str))}
    # Store ks-values for each row for KS test
    ks_values = {row: [] for row in range(len(cobra_reactions_str))}
    # Store fold-change-values for each row for KS test
    fold_change_values = {row: [] for row in range(len(cobra_reactions_str))}

    for row in range(len(cobra_reactions_str)):
        for i, j in selected_comparisons:

            fold_change = np.absolute((np.mean(conditions[i][row]) - np.mean(conditions[j][row])) / (np.mean(conditions[j][row]) + 1e-8))
            # fold_change = (np.mean(conditions[i][row]) - np.mean(conditions[j][row]) ) / (np.mean(conditions[j][row]) + 1e-8)
            fold_change_values[row].append(fold_change)

            if (np.std(conditions[i][row]) > std_cutoff) or (np.std(conditions[j][row]) > std_cutoff):

                ks, p = stats.ks_2samp(conditions[i][row], conditions[j][row], alternative='two-sided')
                p_values[row].append(p)
                ks_values[row].append(ks)

            else:
                p_values[row].append(1)
                ks_values[row].append(1)

    p_values_copy = list(p_values.values())
    flat_p_values = np.array(p_values_copy).flatten()
    # Apply FDR correction
    _, corrected_p_values, _, _ = multipletests(flat_p_values, method='fdr_bh')
    # Reshape corrected p-values back to the original matrix shape
    corrected_p_values = corrected_p_values.reshape(np.array(p_values_copy).shape)

    fold_change_values = np.array(list(fold_change_values.values()))

    significant_diff_indices = np.where(
        np.logical_and(
            corrected_p_values < p_value_cutoff,
            np.abs(fold_change_values) > fold_change_cutoff
        ))[0]

    significant_diff_reactions = [cobra_reactions_str[i] for i in significant_diff_indices]

    not_significant_diff_reactions = list(set(cobra_reactions_str) - set(significant_diff_reactions))

    pval_dict = {cobra_reactions_str[i]: corrected_p_values[i][0] for i in range(len(cobra_reactions_str))}
    fold_change_dict = {cobra_reactions_str[i]: fold_change_values[i][0] for i in range(len(cobra_reactions_str))}

    return significant_diff_reactions, not_significant_diff_reactions, pval_dict, fold_change_dict


def plot_hists(samples, rxn_id, type="dingo", dingo_model=None, description=None, ax=None):
    if type == "dingo":
        samples_df = pd.DataFrame(samples, index=dingo_model.reactions)
        data = samples_df.loc[rxn_id]
    elif type == "cobra":
        data = samples[rxn_id]
    else:
        raise ValueError("Type must be either 'dingo' or 'cobra'.")

    sns.histplot(data, bins=20, kde=False, ax=ax)
    ax.set_xlabel(f"{rxn_id} flux (mmol/gDW/h)")
    ax.set_ylabel("Frequency")
    ax.set_title(f"Histogram of {rxn_id} fluxes {description if description else ''}")


def fba_insight(model, carbon_sources = [], eps = 1e-6):
    """
    For each carbon source, get the reactions that produce or consume it, i.e. flux values that are not zero.
    """
    model.solver = "gurobi"
    sol          = model.optimize()
    sum          = model.summary()

    high_uptakes    = sum.uptake_flux.sort_values("flux", ascending=False).head(10)
    high_secretions = sum.secretion_flux[sum.secretion_flux["flux"] != 0].sort_values("flux")

    if not carbon_sources:
        carbon_source_prefixes = [
            "glc_", "ac_", "lac_", "eth_", "fuc_", "gal_", "fru_",
            "glycerol_", "succ_", "pyruvate_"
        ]
        carbon_sources = set()
        for met in model.metabolites:
            for pre in carbon_source_prefixes:
                if pre in met.id and met.compartment in {"c", "c0", "cytosol"}:
                    carbon_sources.add(met.id)
        carbon_sources = list(carbon_sources)

    cs_rxns = {}
    for csource in carbon_sources:
        met = model.metabolites.get_by_id(csource)
        cs_rxns[csource] = {}
        for rxn in met.reactions:
            v = sol.fluxes[rxn.id]
            # Skip reactions with near-zero flux
            if abs(v) < eps:
                continue
            # Determine if the metabolite is being produced or consumed
            if v > 0:
                if met in rxn.reactants:
                    direction = "consuming"
                else:
                    direction = "producing"
            else:  # v < 0
                if met in rxn.reactants:
                    direction = "producing"
                else:
                    direction = "consuming"

            # # Print with metabolite and reaction information
            # print(
            #     f"{rxn.id} {direction} {met.name} ~ {rxn.build_reaction_string()} || "
            #     f"{rxn.build_reaction_string(use_metabolite_names=True)} *** v = {v:.6f}"
            # )
            cs_rxns[csource][rxn.id] = {
                "direction": direction,
                "metabolite": met.name,
                "reaction": rxn.build_reaction_string(use_metabolite_names=True),
                "flux": f"{v:.2f}"
            }

    return (
        sol,
        high_uptakes,
        high_secretions,
        cs_rxns
    )


# Significantly different fluxes on a dataframe
def different_fluxes_over_sign(samples, rxn_id, reactions_list=None):
    """
    Perform Mann-Whitney U test on each column of the DataFrame to find statistically significant
    """
    if isinstance(samples, np.ndarray):
        if reactions_list is None:
            raise ValueError("If samples is a numpy array, reactions_list must be provided.")
        # Convert numpy array to DataFrame using model reactions as index
        samples = pd.DataFrame(samples, index = reactions_list).T

    # Filter rows where 'col_name' is positive and negative
    positive_rows = samples[samples[rxn_id] > 0]
    negative_rows = samples[samples[rxn_id] < 0]

    print(f"Number of samples with positive flux: {len(positive_rows)}, "
          f"Number of samples with negative flux: {len(negative_rows)}"
          )

    # Initialize lists to store p-values and tested column names
    p_values = []
    tested_columns = []

    # Loop through all columns except 'rxn_id'
    for col in samples.columns:
        if col != rxn_id:
            if np.abs(positive_rows[col].mean()) <= 1e-3 and np.abs(negative_rows[col].mean()) <= 1e-3:
                continue  # skip near-zero flux reactions
            # Perform Mann-Whitney U test
            u_stat, p_value = stats.mannwhitneyu(positive_rows[col], negative_rows[col])
            p_values.append(p_value)
            tested_columns.append(col)

    # Correct for multiple testing
    _, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')

    # Map tested columns to their corrected p-values
    significant_columns = {
        tested_columns[i]: p_values_corrected[i] for i in range(len(p_values_corrected)) if p_values_corrected[i] < 0.05
    }

    print("Statistically significant columns (after FDR correction):")
    print(significant_columns)

    return significant_columns, positive_rows, negative_rows


def count_ios(samples, model: cobra.Model, rxn_id: cobra.Reaction.id = None):
    """
                    0    2    3 
    EX_26dap_M(e)   -    -    +
    EX_ac(e)        +    +    +
    EX_acald(e)     -    +    -

    """
    media_comps = list(model.medium.keys())

    if isinstance(samples, np.ndarray):
        samples = pd.DataFrame(samples, index=[rxn.id for rxn in model.reactions]).T

    columns = []  # store columns as Series
    col_names = []

    for index, sample in samples.iterrows():
        io = sample[media_comps]
        col_values = [
            "+" if v >= 0.01 else "-" if v <= -0.01 else np.nan
            for v in io
        ]

        if not rxn_id or sample[rxn_id] > 0:
            columns.append(pd.Series(col_values, index=media_comps))
            col_names.append(index)

    result_df = pd.concat(columns, axis=1)
    result_df.columns = col_names

    return result_df.T.drop_duplicates().T
