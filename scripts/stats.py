import re
import umap
import cobra
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats
from cobra import Model
from pathlib import Path
from typing import Tuple, List
from numpy.typing import NDArray
from statsmodels.stats.multitest import multipletests

from scipy.stats import spearmanr
import plotly.express as px
import plotly.graph_objects as go
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier


# Apply PCA analysis on a samples dataset
def pca_samples(samples, n_components=10, plots=False, outfile=None):
    """
    samples (pd.DataFrame) -- a df with samples as rows and reactions as column names
    """

    # Perform PCA
    X = StandardScaler().fit_transform(samples)
    pca        = PCA(n_components=n_components)  # You can change the number of components if needed
    pca_result = pca.fit_transform(X)

    # Explained variance ratio
    explained_variance = pca.explained_variance_ratio_
    print(f"Explained Variance Ratio: {explained_variance}")

    if plots:

        if outfile is None:
            outfile = Path() / "pca"
        else:
            outfile = Path(outfile)

        fig, axes = plt.subplots(1, 2, figsize=(12, 6))

        # First subplot: PCA scatter plot
        axes[0].scatter(pca_result[:, 0], pca_result[:, 1])
        axes[0].set_xlabel('Principal Component 1')
        axes[0].set_ylabel('Principal Component 2')
        axes[0].set_xlabel(f'PC1 ({explained_variance[0]*100:.2f}% var)')
        axes[0].set_ylabel(f'PC2 ({explained_variance[1]*100:.2f}% var)')
        axes[0].set_title('PCA: First Two Principal Components')
        axes[0].grid(True)

        # Second subplot: cumulative explained variance
        cumulative_variance = np.cumsum(explained_variance)
        axes[1].plot(range(1, len(cumulative_variance) + 1), cumulative_variance, marker='o')
        axes[1].set_xlabel('Number of Components')
        axes[1].set_ylabel('Cumulative Explained Variance')
        axes[1].set_title('Cumulative Explained Variance by PCA Components')
        axes[1].grid(True)

        # Save the combined figure
        plt.tight_layout()
        plt.savefig(outfile.with_suffix(".eps"), format="eps", dpi=600)
        plt.savefig(outfile.with_suffix(".png"), format="png", dpi=600)

    # --------
    # Get reactions that contribute the most in each PC
    model                = PCA(n_components=2).fit(samples)
    n_pcs                = model.components_.shape[0]
    most_important       = [np.abs(model.components_[i]).argmax() for i in range(n_pcs)]
    most_important_names = [samples.columns[most_important[i]] for i in range(n_pcs)]

    dic = {'PC{}'.format(i): most_important_names[i] for i in range(n_pcs)}

    return dic, model


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

    if ax is None:
        fig, ax = plt.subplots()

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
def different_fluxes_over_sign(samples, rxn_id, model, reactions_list=None):
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

    rows = []
    for rxn_id in significant_columns:
        pos_mean = positive_rows[rxn_id].mean()
        neg_mean = negative_rows[rxn_id].mean()
        if pos_mean * neg_mean < 0:
            rxn = model.reactions.get_by_id(rxn_id)
            rows.append({
                "Reaction ID": rxn.id,
                "Reaction name": rxn.name,
                "Reaction": rxn.build_reaction_string(use_metabolite_names=True),
                "Reactants": [met.id for met in rxn.reactants],
                "Products": [met.id for met in rxn.products],
                "Positive mean": pos_mean,
                "Negative mean": neg_mean
            })

    df_sign_flip = pd.DataFrame(rows)

    return significant_columns, positive_rows, negative_rows, df_sign_flip


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


def tsne(df, perplexity=50, kmeans_k = 2, eps=2, min_samples=5, plot = False, outfile = None):
    """
    t-SNE is a nonlinear dimensionality reduction technique.
    It tries to preserve local neighborhoods: points that are close in high-dimensional space should remain close in 2D/3D.
    """

    # Suppose df is your data matrix (rows = samples, columns = features)
    # Example: df = pd.read_csv("your_data.csv")

    # 1. Standardize features (important for distance-based methods like t-SNE)
    X = StandardScaler().fit_transform(df)

    # 2. (Optional but common) Reduce dimensionality first with PCA
    #    This speeds up t-SNE and denoises
    X_pca = PCA(n_components=50).fit_transform(X)

    # 3. Apply t-SNE
    tsne = TSNE(
        n_components  = 2,    # 2D output
        perplexity    = perplexity,   # typical values: 5–50, tune for your data
        learning_rate = 200,  # also tunable
        max_iter      = 1000,  # Maximum number of iterations for the optimization. Should be at least 250.
        random_state  = 42    # reproducibility
    )
    X_tsne = tsne.fit_transform(X_pca)

    # 4. Put results in a DataFrame for plotting
    tsne_df = pd.DataFrame(X_tsne, columns=["tSNE1", "tSNE2"])

    # Choose a clustering method
    # 1. K-means (requires choosing n_clusters)
    kmeans = KMeans(
        n_clusters   = kmeans_k,
        random_state = 42
    )
    labels_kmeans = kmeans.fit_predict(X_tsne)

    # 2. DBSCAN (density-based, no need to choose k, but need eps & min_samples)
    dbscan = DBSCAN(
        eps         = eps,
        min_samples = min_samples
    )
    labels_dbscan = dbscan.fit_predict(X_tsne)

    if plot:

        if outfile is None:
            outfile = Path() / "dbscan_kmeans_tsne"
        else:
            outfile = Path(outfile)

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # Simple scatter plot
        axes[0].scatter(tsne_df["tSNE1"], tsne_df["tSNE2"], s=30, alpha=0.7)
        axes[0].set_xlabel("t-SNE 1")
        axes[0].set_ylabel("t-SNE 2")
        axes[0].set_title("t-SNE Scatter")

        # Plot K-means result
        axes[1].scatter(X_tsne[:,0], X_tsne[:,1], c=labels_kmeans, cmap='tab10')
        axes[1].set_title('t-SNE + KMeans')

        # Plot DBSCAN result
        axes[2].scatter(X_tsne[:,0], X_tsne[:,1], c=labels_dbscan, cmap='tab10')
        axes[2].set_title('t-SNE + DBSCAN')

        plt.tight_layout()
        plt.show()
        plt.savefig(outfile.with_suffix(".png"), dpi = 300, format="png")

    return labels_kmeans, labels_dbscan


def tsne_dbscan_grid(X, perplexities=[20, 30, 50], min_samples=5, random_state=42):
    """

    Note:

        When to skip PCA before t-SNE:

            If you have moderate dimensions (say ≤50) and the features are already informative.

            If your data is already denoised or you specifically want t-SNE to see the full space.
    """

    X_pca = PCA(n_components=30, random_state=42).fit_transform(X)

    results = {}
    for perp in perplexities:

        # Step 1 – t-SNE projection
        X_tsne = TSNE(n_components=2, perplexity=perp, random_state=random_state).fit_transform(X_pca)

        # Step 2 – k-distance for eps suggestion
        neigh = NearestNeighbors(n_neighbors=min_samples)
        neigh.fit(X_tsne)
        distances, _ = neigh.kneighbors(X_tsne)
        k_distances = np.sort(distances[:, -1])  # distance to k-th neighbor

        # Plot elbow
        plt.figure()
        plt.plot(k_distances)
        plt.title(f"k-distance plot (perplexity={perp})")
        plt.ylabel(f"{min_samples}th NN distance")
        plt.xlabel("Points sorted by distance")
        plt.show()

        # Step 3 – auto-pick eps as elbow
        # Simple heuristic: 90th percentile distance (adjustable)
        eps_guess = np.percentile(k_distances, 90)

        # Step 4 – DBSCAN
        db = DBSCAN(eps=eps_guess, min_samples=min_samples)
        labels = db.fit_predict(X_tsne)

        results[perp] = {
            'tsne'  : X_tsne,
            'eps'   : eps_guess,
            'labels': labels
        }

        print(f"Perplexity={perp}: eps≈{eps_guess:.4f}, clusters={len(set(labels)) - (1 if -1 in labels else 0)}")

    return results


def rxns_cluster_contribution(samples, labels):
    """
    Apply Random Forest tests to get feature contribution on each cluster

    Note:


    samples (pd.DataFraeme)
    labels (np.array)
    """

    df            = samples.copy()
    df["cluster"] = labels

    feature_names = samples.columns

    rf = RandomForestClassifier()
    rf.fit(df[feature_names], df["cluster"])
    importances = pd.Series(rf.feature_importances_, index=feature_names)

    return importances.sort_values(ascending=False)


def is_normal(rxn_distribution):
    """
    Uses the D'Agostino and Pearson's test for normality.
    Returns True if the data is normally distributed, False otherwise.
    """
    from scipy.stats import normaltest
    stat, p_value = normaltest(rxn_distribution)
    print("p =", p_value)
    if p_value > 0.05:
        print("Data appears normally distributed.")
    else:
        print("Data does NOT appear normally distributed.")
    return p_value > 0.05


def check_samples_range(samples, threshold=0.01):

    breaking_constraints = []
    for col in samples:
        if samples[col].min() < -1000 or samples[col].max() > 1000:
            breaking_constraints.append(col)

    non_fixed_fluxes = [
        col for col in samples if abs(samples[col].min() - samples[col].max()) > threshold
    ]
    print(
        f"Number of reactions whose bounds failed: {len(breaking_constraints)}\n"
        f"Number of reaction with a non-fixed flux: {len(non_fixed_fluxes)}"
    )

    return breaking_constraints, non_fixed_fluxes

    "cpd00530_c0",
    "cpd00977_c0",
    "cpd01775_c0"


# Based on the `correlated_reactions()` from:
# github.com/SotirisTouliopoulos/dingo-stats/blob/master/src/correlations_utils.py
def correlated_reactions(
    samples           : pd.DataFrame,
    reactions         : List = [],
    subsystem         : str = None,
    linear_coeff      : str = "pearson",
    linear_corr_cutoff: float = 0.30,
    indicator_cutoff  : float = 1.2,
    lower_triangle    : bool = True,
    label_font_size   : int = 10,
    outfile           : str = None,
    width             : int = 900,
    height            : int = 900
) -> Tuple:
    """
    Correlation matrix for reaction subsystem and according heatmap
    """

    if indicator_cutoff < 1:
        raise Exception("Indicator cutoff must be at least equal to 1")

    samples = samples[reactions].T

    # compute correlation matrix
    if linear_coeff == "pearson":
        corr_matrix = np.corrcoef(samples, rowvar=True)
    elif linear_coeff == "spearman":
        corr_matrix, _ = spearmanr(samples, axis=1)
    else:
        raise Exception("Input value to linear_coeff parameter is not valid. Choose between pearson or spearman")

    # fill diagonal with 1
    np.fill_diagonal(corr_matrix, 1)

    # replace not assigned values with 0
    corr_matrix[np.isnan(corr_matrix)] = 0

    # keep only the lower triangle (with the diagonal) to reduce computational time
    corr_matrix[np.triu_indices(corr_matrix.shape[0], 1)] = np.nan

    # create a copy of correlation matrix to replace/filter values
    correlation_matrix = corr_matrix.copy()

    # find indices of correlation matrix where correlation does not occur
    no_corr_indices = np.argwhere(
        (correlation_matrix < linear_corr_cutoff) & (correlation_matrix > -linear_corr_cutoff)
    )
    # find indices of correlation matrix where correlation does occur
    corr_indices = np.argwhere(
        (correlation_matrix > linear_corr_cutoff) | (correlation_matrix < -linear_corr_cutoff)
    )

    # replace values from the correlation matrix that do not overcome the pearson cutoff with 0
    for i in range(0, no_corr_indices.shape[1]):
        index1 = no_corr_indices[i][0]
        index2 = no_corr_indices[i][1]

        if index1 == index2:
            continue
        else:
            correlation_matrix[index1, index2] = 0
            correlation_matrix[index2, index1] = 0

    correlations_dictionary = {}

    for i in range(0, corr_indices.shape[1]):
        index1 = corr_indices[i][0]
        index2 = corr_indices[i][1]

        if index1 == index2:
            continue

        else:
            reaction1 = reactions[index1]
            reaction2 = reactions[index2]

            pearson = correlation_matrix[index1, index2]

            if pearson > 0:
                classification =  "positive"
            elif pearson < 0:
                classification = "negative"

            correlations_dictionary[reaction1 + "~" + reaction2] = {
                'pearson': pearson,
                'jensenshannon': 0,
                'indicator': 0,
                'classification': classification
            }

    # if user does not want to calculate non-linear correlations
    correlation_matrix[np.triu_indices(correlation_matrix.shape[0], 1)] = np.nan
    np.fill_diagonal(correlation_matrix, 1)

    if not lower_triangle:
        # fill the upper triangle and return a square correlation matrix
        correlation_matrix = np.tril(correlation_matrix) + np.tril(correlation_matrix, -1).T

    # Ensure data is clipped between -1 and +1
    correlation_matrix = np.clip(correlation_matrix, -1, 1)

    sns_colormap = [
        [0.0, '#d73027'],   # red for -1
        [0.5, '#f7f7f7'],   # white for 0
        [1.0, '#4575b4']    # blue for +1
    ]

    fig = px.imshow(
        correlation_matrix,
        color_continuous_scale = sns_colormap,
        zmin                   = -1,
        zmax                   = 1,
        x                      = reactions,
        y                      = reactions,
        origin                 = "upper"
    )

    fig.update_layout(
        xaxis        = dict(tickfont=dict(size=label_font_size)),
        yaxis        = dict(tickfont=dict(size=label_font_size)),
        # x centers the title
        title        = dict(text=f"Correlation Matrix of Reactions: {subsystem}", font=dict(size=20), x=0.5),
        width        = width,
        height       = height,
        plot_bgcolor = "rgba(0,0,0,0)"
    )

    fig.update_traces(xgap=1, ygap=1, hoverongaps=False)
    fig.show()

    safe_subsystem = re.sub(r'_+', "_", re.sub(r'[\\/:"*?<>| ,]+', "_", subsystem))
    if outfile:
        outfile = Path(outfile)
        if not outfile.is_dir():
            dir_path  = outfile.parent
            file_stem = outfile.stem
            outfile = dir_path / f"{file_stem}_cor_heatmap_{safe_subsystem}"
        else:
            outfile = outfile / f"cor_heatmap_{safe_subsystem}"
    else:
        outfile = Path()

    fig.write_image(outfile.with_suffix(".png"), scale=5, width=1200, height=800)
    fig.write_image(outfile.with_suffix(".svg"))

    return correlation_matrix, correlations_dictionary
