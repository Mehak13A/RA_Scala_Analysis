"""
analysis.py

Which genes track the shape of inhibitory neurons in mouse motor cortex?
A single script that runs every step of the analysis, from the three public data files to three result
tables and one figure.

Problem
    Scala et al. (2021) measured gene expression, electrical behaviour and three-dimensional shape in
    the same neurons of mouse motor cortex. Their analysis linked expression to electrical behaviour
    in detail; morphology received less attention. It is not known which genes change expression as a
    neuron's shape changes, or what those genes might do. A naive answer is misleading: neurons come in
    families (Pvalb, Sst, Vip, Lamp5), each with its own gene programme and its own shape, so a gene that
    correlates with shape across all cells may merely be a family marker.

Task
    Identify genes whose expression corresponds with morphological features of inhibitory neurons once
    cell family and position in the cortex are accounted for; interpret their possible function; and
    provide statistical support for the method and the confidence in the result.

Data inputs (from https://github.com/berenslab/mini-atlas, folder data/)
    m1_patchseq_meta_data.csv     : one row per sequenced cell; tab separated
    m1_patchseq_morph_features.csv: one row per reconstructed cell; 63 shape measurements plus depth and class labels
    m1_patchseq_exon_counts.csv.gz: exon read counts; genes as rows, cells as columns

Methodology (each step is a numbered section of the code below)
    Step 1  Select cells. Inhibitory neurons with a reconstructed shape, a reliable family label,
            complete shape measurements and adequate sequencing depth. Every filter is recorded with
            its cell count.
    Step 2  Summarise shape. The 49 shape measurements are standardised and reduced by PCA to five
            shape axes, each named from the measurements that load on it most strongly.
    Step 3  Prepare expression. Counts per million, log2(x + 1), genes detected in too few cells
            removed, and 1,500 highly variable genes selected with the authors' own rule.
    Step 4  Baseline (the authors' approach reproduced). A sparse linear model from genes to shape
            axes with no covariates. Its selected genes are flagged as family markers or not.
    Step 5  Adjusted analysis (the contribution of this work). For every gene and shape axis, a linear
            model with family, soma depth and sequencing depth as covariates, giving an effect size and
            a p-value per gene; false discovery rate controlled with Storey q-values. A multivariate
            comparison of three models (covariates only; genes only; covariates then genes) by cross
            validated R squared shows how much shape genes explain beyond identity and position.
    Step 6  Statistical support. A permutation null for the cross validated R squared (shape rows
            shuffled within family); bootstrap stability of gene selection; and a sensitivity analysis
            on cells with high family label confidence.
    Step 7  Function. Transcription factor activity per cell scored from regulons (decoupler, CollecTRI)
            and tested against the shape axes with the same adjusted model; gene set enrichment of the
            robust genes per axis (g:Profiler). Both need an internet connection and are skipped with a
            message if unavailable.

How to run (from the project folder; config.py in the same folder)
    python analysis.py

Outputs (written to the output folder named in config.py)
    results_gene_associations.csv   one row per gene and shape axis: adjusted effect, p value, q value,
                                    bootstrap stability, whether selected by the baseline model,
                                    whether a family marker
    results_model_summary.csv       per shape axis: cross validated R squared of the three models,
                                    permutation p value, number of genes at q below 0.05, estimated
                                    share of null tests (pi0), sensitivity to label confidence
    results_function.csv            transcription factors and enriched gene sets associated with each
                                    shape axis (empty if step 7 could not run)
    figure_main.png                 six panels summarising the whole analysis
    analysis_log.txt                everything printed during the run, including the cell funnel

References
    1. Scala F, Kobak D, Bernabucci M, et al. Phenotypic variation of transcriptomic cell types in
       mouse motor cortex. Nature 598, 144 to 150 (2021). https://doi.org/10.1038/s41586-020-2907-3
    2. Kobak D, Bernaerts Y, Weis MA, et al. Sparse reduced-rank regression for exploratory visualisation
       of paired multivariate data. J R Stat Soc C 70, 980 to 1000 (2021).
    3. Storey JD, Tibshirani R. Statistical significance for genomewide studies. Proc Natl Acad Sci USA
       100, 9440 to 9445 (2003). Meinshausen N, Buhlmann P. Stability selection. J R Stat Soc B 72, 417
       to 473 (2010).
    4. Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological
       activities from omics data. Bioinform Adv 2, vbac016 (2022).
    5. Muller-Dott S, et al. Expanding the coverage of regulons from high confidence prior knowledge
       for accurate estimation of transcription factor activities. Nucleic Acids Res 51, 10934 to 10949 (2023).
    6. Pedregosa F, et al. Scikit-learn: machine learning in Python. J Mach Learn Res 12, 2825 to 2830 (2011).
"""

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, MultiTaskElasticNet, MultiTaskElasticNetCV
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

warnings.filterwarnings("ignore")
random_state = np.random.RandomState(config.RANDOM_SEED)

class Logger:
    """Print to the screen and to analysis_log.txt at the same time."""

    def __init__(self, path):
        self.file = open(path, "w")

    def __call__(self, text=""):
        print(text)
        self.file.write(str(text) + "\n")
        self.file.flush()


log = Logger(config.OUTPUT_FOLDER / "analysis_log.txt")

# Step 1  Select cells

METADATA_COLUMNS = [
    "RNA family", "RNA type", "RNA type confidence", "Exclusion reasons",
    "Exon reads", "Number of genes detected", "Mouse", "Sequencing batch",
]

def step1_select_cells():
    """
    Join the metadata onto the morphology table, apply the cell filters in a fixed order, and
    return (cell description, shape measurements, exon counts).
    """
    log("STEP 1  Selecting cells")
    metadata = pd.read_csv(config.METADATA_FILE, sep="\t").set_index("Cell")
    morphology = pd.read_csv(config.MORPHOLOGY_FILE, index_col=0)
    funnel = []

    # "keep" function records the count and reason, prints it, and returns the table unchanged
    def keep(table, reason):
        funnel.append((reason, len(table)))
        log(f"  {len(table):5d}  {reason}")
        return table

    cells = keep(morphology, "Reconstructed cells in the morphology table")
    # attach the needed metadata columns to every reconstructed cell using "join"
    cells = cells.join(metadata[METADATA_COLUMNS], how="left")
    cells = keep(cells[cells["cell class"] == config.CELL_CLASS_KEPT], "Inhibitory cells")
    cells = keep(cells[~cells.index.isin(config.ERRATUM_CELLS)], "Reconstruction flagged in the repository errata removed")
    cells = keep(cells[cells["Exclusion reasons"].isna()], "No transcriptomic exclusion reason recorded by the authors")
    cells = keep(cells[cells["RNA family"].isin(config.FAMILIES_KEPT)], f"Family among {config.FAMILIES_KEPT}")
    family_counts = cells["RNA family"].value_counts()
    cells = keep(cells[cells["RNA family"].isin(family_counts[family_counts >= config.MINIMUM_CELLS_PER_FAMILY].index)],
                 f"Family has at least {config.MINIMUM_CELLS_PER_FAMILY} cells")
    cells = keep(cells[cells["normalized depth"].notna()], "Soma depth available")

    # Shape measurements: drop excitatory only columns, position columns and labels.
    drop_columns = [c for c in cells.columns if c.startswith(config.EXCITATORY_ONLY_PREFIX)]
    drop_columns += config.NON_SHAPE_COLUMNS + METADATA_COLUMNS
    shape = cells.drop(columns=drop_columns)
    incomplete = shape.index[shape.isna().any(axis=1)]
    shape = shape.drop(index=incomplete)
    cells = keep(cells.drop(index=incomplete), "Complete set of shape measurements")

    # Exon counts, restricted to these cells and transposed to cells by genes.
    # The count file stores genes as rows and cells as columns, the opposite of every other table here,
    # so it is transposed after the cells are selected.
    counts = pd.read_csv(config.EXON_COUNTS_FILE, index_col=0)
    counts = counts.loc[:, counts.columns.isin(cells.index)].T.loc[cells.index]
    total_reads = counts.sum(axis=1)
    shallow = total_reads.index[total_reads < config.MINIMUM_EXON_READS_PER_CELL]
    counts, shape, cells = counts.drop(index=shallow), shape.drop(index=shallow), cells.drop(index=shallow)
    cells = keep(cells, f"At least {config.MINIMUM_EXON_READS_PER_CELL:,} exon reads in the count matrix")

    description = pd.DataFrame({
        "family": cells["RNA family"],
        "transcriptomic_type": cells["RNA type"],
        "label_confidence": cells["RNA type confidence"],
        "soma_depth": cells["normalized depth"],
        "log10_reads": np.log10(total_reads.loc[cells.index]),
        "genes_detected": cells["Number of genes detected"],
        "mouse": cells["Mouse"],
    })
    log(f"  Kept: {len(description)} neurons, {shape.shape[1]} shape measurements, {counts.shape[1]} genes")
    log("  Cells per family: " + ", ".join(f"{f} {n}" for f, n in description["family"].value_counts().items()))
    return description, shape, counts, funnel

# Step 2  Summarise shape as five axes
# It turns 49 measurements into 5 shape axes and names them.

def step2_shape_axes(shape, description):
    """Standardise, PCA, keep the configured number of axes, name them, check depth."""
    log("\nSTEP 2  Summarising shape")
    standardised = StandardScaler().fit_transform(shape.values) # mean=0, standard deviation=1
    pca = PCA(random_state=config.RANDOM_SEED).fit(standardised) # finds the directions of greatest variation
    number_kept = config.NUMBER_OF_SHAPE_AXES
    axis_ids = [f"shape_axis_{k + 1}" for k in range(number_kept)]
    # give each cell’s score on the first five
    scores = pd.DataFrame(pca.transform(standardised)[:, :number_kept], index=shape.index, columns=axis_ids)
    # give the loadings (how much each measurement contributes to each axis)
    loadings = pd.DataFrame(pca.components_[:number_kept].T, index=shape.columns, columns=axis_ids)
    # give the share of variation per axis
    variance = pca.explained_variance_ratio_

    # A principal component is a weighted sum of all 49 measurements; the weights are the loadings.
    # The sign of a loading says whether the measurement rises or falls along the axis.

    axis_info = []
    for k, axis in enumerate(axis_ids):
        strongest = loadings[axis].abs().sort_values(ascending=False).index[:3]
        words = [("higher " if loadings.loc[m, axis] > 0 else "lower ") + config.readable_label(m).lower() for m in strongest]
        depth_r = stats.pearsonr(scores[axis], description["soma_depth"])[0]
        name = config.SHAPE_AXIS_NAMES.get(axis, ", ".join(words[:2]))
        axis_info.append({"axis": axis, "name": name, "variance_explained": variance[k], "strongest_loadings": "; ".join(words), "depth_correlation": depth_r})
        log(f"  {axis} ({name}): {variance[k]:.1%} of variation, depth r = {depth_r:+.2f}")
    log(f"  Together the {number_kept} axes explain {variance[:number_kept].sum():.1%} of shape variation")
    return scores, loadings, pd.DataFrame(axis_info).set_index("axis"), variance

# Step 3  Prepare expression
# Normalise counts and select 1,500 variable genes.

CLEARLY_ON_COUNT = 32 # the authors' threshold for Patch-seq counts
CURVE_SLOPE = -1.5
CURVE_FLOOR = 0.002

def step3_expression(counts):
    """Counts per million, log2(x + 1), detection filter, authors' variable gene rule."""
    log("\nSTEP 3  Preparing expression")
    # counts in how many cells each gene has at least one read; keep genes with 10 or more
    detected_in = (counts > 0).sum(axis=0)
    counts = counts.loc[:, detected_in >= config.MINIMUM_CELLS_EXPRESSING_GENE]
    log(f"  {counts.shape[1]} genes detected in at least {config.MINIMUM_CELLS_EXPRESSING_GENE} cells")

    # counts per million, which removes the effect of sequencing depth
    log_cpm = np.log2(counts.div(counts.sum(axis=1), axis=0) * config.COUNTS_PER_MILLION_SCALE + 1)

    # The authors' rule: a gene is informative if it is strongly on when on and often off.
    # Genes above a shifted exponential curve are selected; the shift is found by binary search to
    # give the wanted number.
    values = counts.values
    on = values >= CLEARLY_ON_COUNT # marks cells where the gene is clearly on
    number_on = on.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        # mean log2 count over those cells
        mean_when_on = np.where(on, np.log2(np.maximum(values, 1)), 0).sum(axis=0) / np.maximum(number_on, 1)
    fraction_off = 1 - number_on / values.shape[0] # share of cells below 32
    eligible = number_on > 0

    def selected_with(offset):
        return eligible & (fraction_off > np.exp(CURVE_SLOPE * mean_when_on + offset) + CURVE_FLOOR)

    # The offset is found by binary search until exactly 1,500 genes are selected.
    # Raising the offset lifts the curve and selects fewer genes; lowering it selects more.
    # The search halves the interval until the count is exactly the number wanted.
    low, high = -20.0, 20.0
    for _ in range(100):
        middle = (low + high) / 2
        if selected_with(middle).sum() > config.NUMBER_OF_VARIABLE_GENES:
            low = middle
        else:
            high = middle
    selected = counts.columns[selected_with(high)]
    log(f"  {len(selected)} highly variable genes selected with the authors' rule")

    expression = log_cpm[selected]
    standardised = (expression - expression.mean()) / expression.std(ddof=0)
    return log_cpm, expression, standardised

# Shared helpers for the models

# a short hand-written list of well-known family marker genes, used only to flag genes in the output
FAMILY_MARKERS = ["Pvalb", "Sst", "Vip", "Lamp5", "Sncg", "Gad1", "Gad2", "Slc6a1",
                  "Npy", "Calb1", "Calb2", "Reln", "Tac1", "Htr3a", "Cck", "Crh", "Nos1", "Chodl"]

def covariate_matrix(description): # builds the table of covariates
    """Family (as indicator columns), soma depth and sequencing depth."""
    family_dummies = pd.get_dummies(description["family"], prefix="family", drop_first=True).astype(float)
    covariates = pd.concat([family_dummies, description[["soma_depth", "log10_reads"]]], axis=1)
    return covariates

def cross_validated_predictions(y, predict_function):
    """K fold cross validated predictions for every row of y."""
    # split the cells into 10 groups; for each group, the model is fitted on the other 9 and predicts this one
    folds = KFold(config.CROSS_VALIDATION_FOLDS, shuffle=True, random_state=config.RANDOM_SEED)
    predictions = np.zeros_like(y, dtype=float)
    for train, test in folds.split(y):
        predictions[test] = predict_function(train, test)
    return predictions

def cross_validated_r2(y, predict_function):
    """
    Generic K fold cross validation. predict_function(train_index, test_index) returns predictions
    for the test rows. Returns R squared per column of y.
    """
    folds = KFold(config.CROSS_VALIDATION_FOLDS, shuffle=True, random_state=config.RANDOM_SEED)
    predictions = np.zeros_like(y, dtype=float)
    for train, test in folds.split(y):
        predictions[test] = predict_function(train, test)
    residual = ((y - predictions) ** 2).sum(axis=0)
    total = ((y - y.mean(axis=0)) ** 2).sum(axis=0)
    return 1 - residual / total

def storey_q_values(p_values):
    """
    Storey and Tibshirani (2003). pi0, the share of tests with no real effect, is estimated from the
    flat part of the p value histogram; q-values are then Benjamini Hochberg values scaled by pi0.
    """
    p = np.asarray(p_values)
    m = len(p)
    lambdas = np.arange(0.05, 0.96, 0.05)
    pi0_by_lambda = np.array([(p > lam).mean() / (1 - lam) for lam in lambdas])
    # Smooth with a cubic spline and take the value at the largest lambda.
    from scipy.interpolate import UnivariateSpline
    spline = UnivariateSpline(lambdas, pi0_by_lambda, k=3, s=len(lambdas))
    pi0 = float(np.clip(spline(lambdas[-1]), 0.0, 1.0))
    if pi0 == 0:
        pi0 = 1e-6
    order = np.argsort(p)
    ranks = np.arange(1, m + 1)
    q_sorted = pi0 * m * p[order] / ranks
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.minimum(q_sorted, 1.0)
    return q, pi0

def gene_wise_associations(axes, expression_std, covariates): # gene-by-gene adjusted model
    """
    For every gene and axis: slope, p value of the gene term in
        axis ~ gene + family + soma depth + sequencing depth.
    Implemented with the Frisch Waugh Lovell trick: residualise both the axis and every gene on the
    covariates, then the slope of residual on residual equals the gene coefficient of the full model.
    Degrees of freedom account for the covariates.
    """

    # Build the projection matrix I - X(X'X)^-1 X' from the covariates, apply it to the axes and to
    # every gene (this removes what the covariates explain from both), then the slope of residual axis
    # on residual gene equals the gene’s coefficient in the full model.

    X = np.column_stack([np.ones(len(covariates)), covariates.values])
    projection = np.eye(len(X)) - X @ np.linalg.pinv(X)
    axes_residual = projection @ axes.values
    genes_residual = projection @ expression_std.values
    degrees_of_freedom = len(X) - X.shape[1] - 1

    rows = []
    for a, axis in enumerate(axes.columns):
        y = axes_residual[:, a]
        gxx = (genes_residual ** 2).sum(axis=0)
        slope = (genes_residual * y[:, None]).sum(axis=0) / gxx
        residual_variance = ((y[:, None] - genes_residual * slope) ** 2).sum(axis=0) / degrees_of_freedom
        standard_error = np.sqrt(residual_variance / gxx)
        t = slope / standard_error
        p = 2 * stats.t.sf(np.abs(t), degrees_of_freedom)
        rows.append(pd.DataFrame({"gene": expression_std.columns, "axis": axis, "adjusted_effect": slope, "t_statistic": t, "p_value": p}))
    return pd.concat(rows, ignore_index=True)

# Grid of penalty strengths searched by cross validation. The data are standardised, so the useful
# range is known; written explicitly because the argument that generates a grid automatically differs
# between scikit-learn versions.
PENALTY_GRID = np.logspace(-1.5, 1.0, 30)

# Step 4  Baseline: the authors' approach reproduced

def step4_baseline(axes, expression_std):
    """
    Sparse multi output linear model from all genes to all shape axes with no covariates
    (the same model family as the authors' sparse reduced rank regression: linear, penalised, few genes
    kept). Returns the selected genes and their cross validated R squared per axis.
    """
    log("\nSTEP 4  Baseline: sparse regression from genes to shape, no covariates")
    y = axes.values
    # a linear model predicting all five axes at once from all 1,500 genes, with a penalty that shrinks
    # coefficients and sets most to exactly zero (“sparse”)

    # "l1_ratio=0.5" mixes lasso (sparsity) and ridge (stability).

    # "cv" means it tries each penalty in the grid and picks the one with the best cross-validated error
    model = MultiTaskElasticNetCV(l1_ratio=0.5, cv=config.CROSS_VALIDATION_FOLDS, random_state=config.RANDOM_SEED, alphas=PENALTY_GRID, max_iter=5000).fit(expression_std.values, y)
    alpha = model.alpha_
    # selected genes are whose coefficient is non-zero on any axis.
    selected = expression_std.columns[np.abs(model.coef_).sum(axis=0) > 0]
    r2 = cross_validated_r2(y, lambda tr, te: MultiTaskElasticNet(alpha=alpha, l1_ratio=0.5, max_iter=5000).fit(expression_std.values[tr], y[tr]).predict(expression_std.values[te]))
    markers = [g for g in selected if g in FAMILY_MARKERS]
    log(f"  {len(selected)} genes selected; {len(markers)} are known family markers: {markers}")
    log("  Cross validated R squared per axis: " + ", ".join(f"{v:.2f}" for v in r2))
    return list(selected), r2, alpha

# Step 5  Adjusted analysis

def step5_adjusted(axes, expression_std, description, baseline_alpha):
    """
    Gene wise adjusted associations with Storey q values, and the three-way cross validated comparison:
    covariates only; genes only; covariates then genes on the residual.
    """
    # Calls gene_wise_associations with the covariates; applies storey_q_values across all 7,500 p-values;
    # counts significant genes per axis; prints the ten strongest.

    log("\nSTEP 5  Adjusted analysis: family, soma depth and sequencing depth as covariates")
    covariates = covariate_matrix(description)
    associations = gene_wise_associations(axes, expression_std, covariates)
    associations["q_value"], pi0 = storey_q_values(associations["p_value"].values)
    significant = associations[associations["q_value"] < config.FALSE_DISCOVERY_THRESHOLD]
    log(f"  Estimated share of gene by axis tests with no real effect (pi0): {pi0:.2f}")
    log(f"  Gene by axis associations at q < {config.FALSE_DISCOVERY_THRESHOLD}: {len(significant)}")
    for axis in axes.columns:
        log(f"    {axis}: {int((significant['axis'] == axis).sum())} genes")
    strongest = associations.nsmallest(10, "p_value")
    log("  Ten strongest adjusted associations (gene, axis, effect, p, q):")
    for _, row in strongest.iterrows():
        log(f"    {row['gene']:12s} {row['axis']}  effect {row['adjusted_effect']:+.3f}  p = {row['p_value']:.2e}  q = {row['q_value']:.3f}")

    y, C, G = axes.values, covariates.values, expression_std.values

    # The penalty for the gene model on the residuals is tuned on the residuals themselves: they have
    # much less variance than the raw axes, so reusing the baseline penalty would shrink every gene to zero.
    residual_full = y - LinearRegression().fit(C, y).predict(C)
    residual_alpha = MultiTaskElasticNetCV(l1_ratio=0.5, cv=config.CROSS_VALIDATION_FOLDS, random_state=config.RANDOM_SEED, alphas=PENALTY_GRID, max_iter=5000).fit(G, residual_full).alpha_
    log(f"  Penalty chosen for the gene model on residuals: {residual_alpha:.4f} (baseline penalty {baseline_alpha:.4f})")

    def covariates_only(tr, te): # ordinary least squares on family and depth
        return LinearRegression().fit(C[tr], y[tr]).predict(C[te])

    def genes_only(tr, te): # elastic net with the baseline penalty
        return MultiTaskElasticNet(alpha=baseline_alpha, l1_ratio=0.5, max_iter=5000).fit(G[tr], y[tr]).predict(G[te])

    def covariates_then_genes(tr, te):
        # fit covariates, take the residual on the training fold, fit the gene model on the residual,
        # add the two predictions
        base = LinearRegression().fit(C[tr], y[tr])
        residual_train = y[tr] - base.predict(C[tr])
        gene_model = MultiTaskElasticNet(alpha=residual_alpha, l1_ratio=0.5, max_iter=5000).fit(G[tr], residual_train)
        return base.predict(C[te]) + gene_model.predict(G[te])

    predictions_covariates = cross_validated_predictions(y, covariates_only)
    predictions_genes = cross_validated_predictions(y, genes_only)
    r2_covariates = cross_validated_r2(y, covariates_only)
    r2_genes = cross_validated_r2(y, genes_only)
    r2_both = cross_validated_r2(y, covariates_then_genes)

    # How similar are the two sets of predictions? If the genes only model predicts the same values as
    # the family and depth model, then what the genes "know" about shape is family identity and position.
    prediction_agreement = np.array([stats.pearsonr(predictions_genes[:, k], predictions_covariates[:, k])[0] for k in range(y.shape[1])])
    log("  Agreement between genes only and family and depth only predictions (r): " + ", ".join(f"{v:.2f}" for v in prediction_agreement))
    for k, axis in enumerate(axes.columns):
        log(f"  {axis}: R2 covariates only {r2_covariates[k]:.2f}, genes only {r2_genes[k]:.2f}, "
            f"covariates then genes {r2_both[k]:.2f}, gain from genes {r2_both[k] - r2_covariates[k]:+.2f}")
    return associations, pi0, r2_covariates, r2_genes, r2_both, covariates, residual_alpha, prediction_agreement

# Step 5b  Within family continuity at the cell level

def step5b_continuity(axes, expression_std, description, covariates):
    """
    The paper tested continuity by correlating distances between type averages. Here the same question
    is asked at the level of individual cells, within each family: is a pair of cells that is close in
    gene space also close in shape space? Gene distance is one minus the Pearson correlation between
    two cells across the 1,500 genes. Shape distance is the Euclidean distance between the cells on the
    five axes after family, depth and sequencing depth have been removed (so position does not drive
    the result). The correlation between the two distance vectors is tested with a Mantel permutation
    (cell labels of one matrix shuffled).
    """
    log("\nSTEP 5b  Within family continuity: do transcriptomic neighbours have similar shapes?")
    y, C = axes.values, covariates.values
    residual_axes = y - LinearRegression().fit(C, y).predict(C)
    residual_axes = (residual_axes - residual_axes.mean(axis=0)) / residual_axes.std(axis=0)
    results = {}
    for family in config.FAMILIES_KEPT:
        rows = np.where(description["family"].values == family)[0]
        if len(rows) < config.MINIMUM_CELLS_PER_FAMILY:
            continue
        # give the gene correlation between every pair of cells
        gene_corr = np.corrcoef(expression_std.values[rows])
        gene_distance = 1 - gene_corr
        shape_distance = np.sqrt(((residual_axes[rows][:, None, :] - residual_axes[rows][None, :, :]) ** 2).sum(axis=2))
        upper = np.triu_indices(len(rows), k=1) # picks each pair once
        # correlate the two distance vectors (the Mantel statistic)
        observed = stats.spearmanr(gene_distance[upper], shape_distance[upper])[0]
        null = np.empty(config.NUMBER_OF_PERMUTATIONS)
        for i in range(config.NUMBER_OF_PERMUTATIONS):
            order = random_state.permutation(len(rows))
            null[i] = stats.spearmanr(gene_distance[upper], shape_distance[order][:, order][upper])[0]
        p = (np.sum(null >= observed) + 1) / (config.NUMBER_OF_PERMUTATIONS + 1)
        results[family] = {"cells": len(rows), "mantel_r": observed, "null_95th": np.percentile(null, 95), "p": p}
        log(f"  {family}: {len(rows)} cells, Mantel r = {observed:+.3f}, null 95th percentile {np.percentile(null, 95):+.3f}, p = {p:.3f}")
    return pd.DataFrame(results).T

# Step 6  Statistical support

def step6_support(axes, expression_std, description, covariates, associations, residual_alpha):
    """Permutation null for the gain from genes, bootstrap stability, confidence sensitivity."""
    log("\nSTEP 6  Statistical support")
    y, C, G = axes.values, covariates.values, expression_std.values
    base_fit = LinearRegression().fit(C, y)
    residual_axes = y - base_fit.predict(C)
    families = description["family"].values

    # Permutation: shuffle residual shape rows within each family, refit the gene model on the
    # residuals with the chosen penalty, record cross validated R squared of the gene part.
    # The observed gain is compared to this null distribution.
    def gene_r2_on(residual):
        return cross_validated_r2(residual, lambda tr, te: MultiTaskElasticNet(alpha=residual_alpha, l1_ratio=0.5, max_iter=3000).fit(G[tr], residual[tr]).predict(G[te]))

    observed_gain = gene_r2_on(residual_axes)
    null_gains = np.zeros((config.NUMBER_OF_PERMUTATIONS, axes.shape[1]))
    for i in range(config.NUMBER_OF_PERMUTATIONS):
        shuffled = residual_axes.copy()
        for family in np.unique(families):
            rows = np.where(families == family)[0]
            shuffled[rows] = shuffled[random_state.permutation(rows)]
        null_gains[i] = gene_r2_on(shuffled)
        if (i + 1) % 50 == 0:
            log(f"    permutation {i + 1} of {config.NUMBER_OF_PERMUTATIONS}")
    permutation_p = ((null_gains >= observed_gain).sum(axis=0) + 1) / (config.NUMBER_OF_PERMUTATIONS + 1)
    for k, axis in enumerate(axes.columns):
        log(f"  {axis}: gene R2 on residual {observed_gain[k]:.3f}, null 95th percentile "
            f"{np.percentile(null_gains[:, k], 95):.3f}, permutation p = {permutation_p[k]:.3f}")

    # Bootstrap stability: resample cells, refit the sparse gene model on the residuals, record which
    # genes are selected. A gene selected in most resamples is robust to which cells happened to be measured.
    selection_counts = pd.Series(0.0, index=expression_std.columns)
    for i in range(config.NUMBER_OF_BOOTSTRAPS):
        sample = random_state.choice(len(y), len(y), replace=True)
        fit = MultiTaskElasticNet(alpha=residual_alpha, l1_ratio=0.5, max_iter=3000).fit(G[sample], residual_axes[sample])
        selection_counts += (np.abs(fit.coef_).sum(axis=0) > 0)
    # stability is the selection frequency per gene; robust genes exceed 70%
    stability = selection_counts / config.NUMBER_OF_BOOTSTRAPS
    robust_genes = stability.index[stability >= config.STABILITY_THRESHOLD].tolist()
    log(f"  Genes selected in at least {config.STABILITY_THRESHOLD:.0%} of {config.NUMBER_OF_BOOTSTRAPS} bootstraps: {len(robust_genes)}")

    # Sensitivity to label confidence: rerun the gene wise model on confident cells only and correlate
    # the effects with the full run.
    sensitivity = {}
    for threshold in config.LABEL_CONFIDENCE_THRESHOLDS:
        confident = description["label_confidence"] >= threshold
        sub_assoc = gene_wise_associations(axes[confident], expression_std[confident], covariates[confident])
        merged = associations.merge(sub_assoc, on=["gene", "axis"], suffixes=("", "_confident"))
        top = merged.nsmallest(100, "p_value")
        r = stats.pearsonr(top["adjusted_effect"], top["adjusted_effect_confident"])[0]
        sensitivity[threshold] = (int(confident.sum()), r)
        log(f"  Label confidence >= {threshold}: {int(confident.sum())} cells; effect sizes of the 100 strongest "
            f"gene by axis associations correlate with the full run at r = {r:.2f}")
    return observed_gain, null_gains, permutation_p, stability, robust_genes, sensitivity

# Step 7  Function (needs internet; skipped gracefully if unavailable)

def step7_function(axes, log_cpm, description, covariates, associations, robust_genes):
    """Transcription factor activity against shape axes, and enrichment of robust genes."""
    log("\nSTEP 7  Function")
    results = []

    # Transcription factor activity with decoupler on the CollecTRI network.
    try:
        import decoupler as dc
        # downloads the regulon network (which TF regulates which genes)
        network = dc.op.collectri(organism="mouse")
        # fits, per cell and per TF, a small linear model asking whether that TF’s targets are
        # collectively up or down; the result is a cells-by-TFs activity table
        activity = dc.mt.ulm(data=log_cpm, net=network, tmin=5)
        tf_scores = activity[0] if isinstance(activity, tuple) else activity
        tf_std = (tf_scores - tf_scores.mean()) / tf_scores.std(ddof=0)
        tf_assoc = gene_wise_associations(axes, tf_std, covariates)
        tf_assoc["q_value"], _ = storey_q_values(tf_assoc["p_value"].values)
        tf_assoc = tf_assoc.rename(columns={"gene": "item"})
        tf_assoc["kind"] = "transcription_factor"
        significant_tf = tf_assoc[tf_assoc["q_value"] < config.FALSE_DISCOVERY_THRESHOLD]
        results.append(significant_tf)
        log(f"  Transcription factors scored: {tf_scores.shape[1]}; associations at q < 0.05: {len(significant_tf)}")
        strongest_tf = tf_assoc.nsmallest(5, "p_value")
        log("  Five strongest transcription factor associations (factor, axis, effect, p, q):")
        for _, row in strongest_tf.iterrows():
            log(f"    {row['item']:12s} {row['axis']}  effect {row['adjusted_effect']:+.3f}  p = {row['p_value']:.2e}  q = {row['q_value']:.3f}")
        # The strongest five are kept in the results table even when none reaches q < 0.05, so the
        # table records what was tested.
        strongest_tf = strongest_tf.copy()
        strongest_tf["kind"] = "transcription_factor_strongest_five"
        results.append(strongest_tf)
    except Exception as error:
        log(f"  Transcription factor step skipped ({type(error).__name__}: {error})")

    # Enrichment of the robust genes, per axis and direction, with g:Profiler.
    # For each axis and direction, take genes that are robust and significant, query g:Profiler against
    # the 1,500 tested genes as background.
    try:
        from gprofiler import GProfiler
        gp = GProfiler(return_dataframe=True)
        background = associations["gene"].unique().tolist()
        for axis in axes.columns:
            for direction, sign in [("increases_with_axis", 1), ("decreases_with_axis", -1)]:
                subset = associations[(associations["axis"] == axis) & (associations["gene"].isin(robust_genes)) & (np.sign(associations["adjusted_effect"]) == sign) & (associations["q_value"] < config.FALSE_DISCOVERY_THRESHOLD)]
                if len(subset) < 5:
                    continue
                enriched = gp.profile(organism="mmusculus", query=subset["gene"].tolist(), background=background, sources=["GO:BP", "KEGG", "REAC"], no_evidences=False)
                if len(enriched):
                    enriched = enriched[["name", "p_value", "intersection_size", "source"]].head(10)
                    enriched = enriched.rename(columns={"name": "item"})
                    enriched["axis"], enriched["kind"] = axis, f"enriched_set_{direction}"
                    results.append(enriched)
        enriched_rows = sum(len(r) for r in results if "kind" in r.columns and r["kind"].astype(str).str.startswith("enriched").any())
        log(f"  Enrichment completed: {enriched_rows} enriched gene sets written "
            f"(enrichment runs only on genes that are both bootstrap stable and at q < 0.05)")
    except Exception as error:
        log(f"  Enrichment step skipped ({type(error).__name__}: {error})")

    return pd.concat(results, ignore_index=True) if results else pd.DataFrame(columns=["item", "axis", "kind"])

# The figure

def draw_figure(description, axis_info, variance, baseline_genes, associations, prediction_agreement,
                r2_covariates, r2_genes, r2_both, observed_gain, null_gains, permutation_p, continuity, path):
    axis_ids = list(axis_info.index)
    short_names = [f"Axis {k + 1}: {axis_info.loc[a, 'name']}" for k, a in enumerate(axis_ids)]
    figure, panels = plt.subplots(2, 3, figsize=(20, 12))
    figure.suptitle("Which genes track the shape of inhibitory neurons? Scala et al. 2021 Patch-seq data, "
                    f"{len(description)} neurons", fontsize=15)

    # (a) cells per family and their depth.
    panel = panels[0, 0]
    families = [f for f in config.FAMILIES_KEPT if f in description["family"].unique()]
    for i, family in enumerate(families):
        depths = description.loc[description["family"] == family, "soma_depth"]
        panel.scatter(np.random.normal(i, 0.07, len(depths)), depths, s=12, alpha=0.6, color=config.FAMILY_COLOURS[family])
    panel.set_xticks(range(len(families)))
    panel.set_xticklabels([f"{f}\nn = {int((description['family'] == f).sum())}" for f in families])
    panel.invert_yaxis()
    panel.set_ylabel(config.readable_label("normalized depth"))
    panel.set_title("a. The cells: four families, and where they sit")

    # (b) the five shape axes.
    panel = panels[0, 1]
    shown = min(12, len(variance))
    panel.bar(np.arange(1, shown + 1), variance[:shown] * 100, color="#5d6d7e")
    panel.plot(np.arange(1, shown + 1), np.cumsum(variance[:shown]) * 100, marker="o", color="#c0392b")
    panel.axvline(len(axis_ids) + 0.5, linestyle=":", color="black")
    panel.text(0.98, 0.04, "\n".join(short_names), transform=panel.transAxes, ha="right", va="bottom", fontsize=9, bbox=dict(boxstyle="round", facecolor="white", edgecolor="#aeb6bf"))
    panel.set_xlabel("Shape axis (principal component)")
    panel.set_ylabel("Share of shape variation explained (%)")
    panel.set_title("b. Shape summarised as five axes")

    # (c) how much shape genes explain, with the permutation null.
    panel = panels[0, 2]
    positions = np.arange(len(axis_ids))
    width = 0.25
    panel.bar(positions - width, r2_covariates, width, color="#aeb6bf", label="Family and depth only")
    panel.bar(positions, r2_genes, width, color="#5dade2", label="Genes only (authors' approach)")
    panel.bar(positions + width, r2_both, width, color="#1a5276", label="Family and depth, then genes")
    for k in range(len(axis_ids)):
        panel.text(positions[k] + width, r2_both[k] + 0.01, f"p = {permutation_p[k]:.2f}", ha="center", fontsize=8)
    panel.set_xticks(positions)
    panel.set_xticklabels([f"Axis {k + 1}" for k in range(len(axis_ids))])
    panel.set_ylabel("Prediction accuracy (cross validated R squared)")
    panel.set_title("c. How much shape genes explain beyond family and depth\n(p from permutation of shape within family)")
    panel.legend(fontsize=8, loc="upper right")

    # (d) what the genes predict is what family and depth predict.
    panel = panels[1, 0]
    positions = np.arange(len(axis_ids))
    panel.bar(positions, prediction_agreement, color="#1a5276")
    for k, value in enumerate(prediction_agreement):
        panel.text(positions[k], value + 0.01, f"{value:.2f}", ha="center", fontsize=9)
    panel.set_ylim(0, 1.05)
    panel.set_xticks(positions)
    panel.set_xticklabels([f"Axis {k + 1}" for k in range(len(axis_ids))])
    panel.set_ylabel("Correlation between the two models' predictions")
    panel.set_title(f"d. The genes only model ({len(baseline_genes)} genes) predicts\nthe same values as family and depth alone")

    # (e) heatmap of the robust genes against the axes (adjusted effect, stars at q < 0.05).
    panel = panels[1, 1]
    top_genes = associations.nsmallest(25, "p_value")["gene"].unique()
    top = associations[associations["gene"].isin(top_genes)]
    if len(top):
        order = (top.groupby("gene")["p_value"].min().sort_values().index)
        matrix = top.pivot(index="gene", columns="axis", values="adjusted_effect").loc[order, axis_ids]
        qmat = top.pivot(index="gene", columns="axis", values="q_value").loc[order, axis_ids]
        limit = np.nanmax(np.abs(matrix.values))
        image = panel.imshow(matrix.values, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                if qmat.values[i, j] < config.FALSE_DISCOVERY_THRESHOLD:
                    panel.text(j, i, "*", ha="center", va="center", fontsize=9)
        panel.set_yticks(range(len(order)))
        panel.set_yticklabels(order, fontsize=8)
        panel.set_xticks(range(len(axis_ids)))
        panel.set_xticklabels([f"Axis {k + 1}" for k in range(len(axis_ids))])
        figure.colorbar(image, ax=panel, fraction=0.04, label="Adjusted effect (standard deviations of axis per SD of gene)")
    min_q = associations["q_value"].min()
    panel.set_title(f"e. Strongest adjusted associations (* q < 0.05; smallest q = {min_q:.2f})")

    # (f) within family continuity at the cell level.
    panel = panels[1, 2]
    families_tested = list(continuity.index)
    positions = np.arange(len(families_tested))
    panel.bar(positions, continuity["mantel_r"], color=[config.FAMILY_COLOURS[f] for f in families_tested])
    panel.scatter(positions, continuity["null_95th"], marker="_", s=600, color="black", label="95th percentile of shuffled null")
    for k, family in enumerate(families_tested):
        panel.text(positions[k], continuity.loc[family, "mantel_r"] + 0.005,
                   f"p = {continuity.loc[family, 'p']:.2f}\nn = {int(continuity.loc[family, 'cells'])}", ha="center", fontsize=8)
    panel.axhline(0, color="black", linewidth=0.8)
    panel.set_ylim(min(0, continuity["mantel_r"].min() - 0.02), max(continuity["mantel_r"].max(), continuity["null_95th"].max()) * 1.6)
    panel.set_xticks(positions)
    panel.set_xticklabels(families_tested)
    panel.set_ylabel("Mantel correlation: gene distance versus shape distance between cells")
    panel.set_title("f. Within a family, are transcriptomic neighbours shape neighbours?")
    panel.legend(fontsize=8, loc="upper right")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    figure.savefig(path, dpi=150)
    plt.close(figure)
    log(f"\nFigure written: {path.name}")

# Main

def main():
    description, shape, counts, funnel = step1_select_cells()
    axes, loadings, axis_info, variance = step2_shape_axes(shape, description)
    log_cpm, expression, expression_std = step3_expression(counts)
    baseline_genes, r2_baseline, alpha = step4_baseline(axes, expression_std)
    associations, pi0, r2_covariates, r2_genes, r2_both, covariates, residual_alpha, prediction_agreement = step5_adjusted(
        axes, expression_std, description, alpha)
    continuity = step5b_continuity(axes, expression_std, description, covariates)
    observed_gain, null_gains, permutation_p, stability, robust_genes, sensitivity = step6_support(
        axes, expression_std, description, covariates, associations, residual_alpha)
    function_table = step7_function(axes, log_cpm, description, covariates, associations, robust_genes)

    # Result table 1: gene associations.
    associations["bootstrap_stability"] = associations["gene"].map(stability)
    associations["selected_by_baseline"] = associations["gene"].isin(baseline_genes)
    associations["known_family_marker"] = associations["gene"].isin(FAMILY_MARKERS)
    associations["axis_name"] = associations["axis"].map(axis_info["name"])
    associations = associations.sort_values(["q_value", "p_value"])
    associations.to_csv(config.OUTPUT_FOLDER / "results_gene_associations.csv", index=False)

    # Result table 2: model summary per axis.
    summary = axis_info.copy()
    summary["r2_covariates_only"] = r2_covariates
    summary["r2_genes_only_baseline"] = r2_genes
    summary["r2_covariates_then_genes"] = r2_both
    summary["agreement_genes_vs_covariate_predictions"] = prediction_agreement
    summary["gene_r2_on_residual"] = observed_gain
    summary["permutation_null_95th_percentile"] = np.percentile(null_gains, 95, axis=0)
    summary["permutation_p_value"] = permutation_p
    summary["genes_at_q_below_threshold"] = [int(((associations["axis"] == a) & (associations["q_value"] < config.FALSE_DISCOVERY_THRESHOLD)).sum()) for a in summary.index]
    summary["pi0_share_of_null_tests"] = pi0
    summary["baseline_genes_selected"] = len(baseline_genes)
    summary["robust_genes"] = len(robust_genes)
    for threshold, (n_cells, r) in sensitivity.items():
        summary[f"effect_correlation_confidence_{threshold}"] = r
        summary[f"cells_confidence_{threshold}"] = n_cells
    for family, row in continuity.iterrows():
        summary[f"continuity_mantel_r_{family}"] = row["mantel_r"]
        summary[f"continuity_p_{family}"] = row["p"]
    summary.to_csv(config.OUTPUT_FOLDER / "results_model_summary.csv")

    # Result table 3: function.
    function_table.to_csv(config.OUTPUT_FOLDER / "results_function.csv", index=False)

    draw_figure(description, axis_info, variance, baseline_genes, associations, prediction_agreement,
                r2_covariates, r2_genes, r2_both, observed_gain, null_gains, permutation_p, continuity,
                config.OUTPUT_FOLDER / "figure_main.png")
    log("\nDone. Outputs: results_gene_associations.csv, results_model_summary.csv, results_function.csv, figure_main.png")

if __name__ == "__main__":
    main()