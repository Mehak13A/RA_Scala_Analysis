"""
config.py
Shared settings for the analysis of gene expression and neuron morphology.

Problem
    Scala et al. (2021) measured gene expression, electrical behaviour and three-dimensional shape
    in the same neurons of mouse primary motor cortex. The published analysis linked gene expression
    to electrical behaviour in detail but treated morphology more briefly. It remains an open question
    which genes change their expression as the shape of a neuron changes, and whether those genes point
    to a biological function.

Task
    Using the publicly released Patch-seq data, identify genes whose expression corresponds with
    morphological features of inhibitory neurons, interpret the possible function of those genes,
    and provide statistical support for the method and the confidence in the results. The analysis
    is restricted to inhibitory neurons because they are the only class for which both axon and
    dendrite were reconstructed.

Role of this file
    Every decision that could reasonably be made differently lives here, so that the reviewer can see
    each choice in one place and change it in one place. The numbered scripts contain the analysis and
    import their settings from this file. Nothing in this file performs any computation.

Data inputs (from https://github.com/berenslab/mini-atlas, folder data/)
    m1_patchseq_meta_data.csv     : one row per sequenced cell; tab separated
    m1_patchseq_morph_features.csv: one row per reconstructed cell; 63 shape measurements plus depth and class labels
    m1_patchseq_exon_counts.csv.gz: exon read counts; genes as rows, cells as columns

How to use
    This file and analysis.py sit together in the project folder, alongside  a "data" folder holding
    the three mini-atlas files and an "outputs" folder where results are written. The two folder paths
    in the section "Folder locations" are found from the location of this file, so the analysis runs
    from any copy of the project folder without editing. To use other folders, the two lines can be
    replaced with explicit paths, as shown in the comment above them. The analysis is then run from the
    project folder with: python analysis.py

Outputs
    None directly. All outputs are written into the folder named in OUTPUT_DIR.

References
    1. Scala F, Kobak D, Bernabucci M, et al. Phenotypic variation of transcriptomic cell types in mouse
       motor cortex. Nature 598, 144 to 150 (2021). https://doi.org/10.1038/s41586-020-2907-3
    2. Berens laboratory. mini-atlas repository: analysis code and preprocessed data for Scala et al.
       2021. https://github.com/berenslab/mini-atlas
"""

from pathlib import Path

# Folder locations
# PROJECT_ROOT is the folder containing this file. DATA_DIR holds the three mini-atlas files and
# OUTPUT_DIR receives every table, figure and log. Both are found relative to this file, so no path
# needs changing on another computer. To point elsewhere, replace the two lines with explicit paths,
# for example on a Mac:
#     DATA_DIR = "/Users/name/ra_scala_interview_2/data"
#     OUTPUT_DIR = "/Users/name/ra_scala_interview_2/outputs"
# or on Windows (the r before the quote keeps the backslashes literal):
#     DATA_DIR = r"C:\Users\name\ra_scala_interview_2\data"
#     OUTPUT_DIR = r"C:\Users\name\ra_scala_interview_2\outputs"
PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data" # inputs
OUTPUT_DIR = PROJECT_ROOT / "outputs" # outputs

DATA_FOLDER = Path(DATA_DIR)
OUTPUT_FOLDER = Path(OUTPUT_DIR)
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

# The three mini-atlas files, located inside the data folder.
METADATA_FILE = DATA_FOLDER / "m1_patchseq_meta_data.csv"
MORPHOLOGY_FILE = DATA_FOLDER / "m1_patchseq_morph_features.csv"
EXON_COUNTS_FILE = DATA_FOLDER / "m1_patchseq_exon_counts.csv.gz"

# Random seed
# A single integer used by every script wherever randomness is involved (bootstraps, permutations,
# cross validation folds), so a rerun gives the same numbers.
RANDOM_SEED = 20261030

# Cell selection
# The morphology table labels each cell "inh" (inhibitory) or "exc" (excitatory).
# Only inhibitory cells are analysed, for the reason given in the Task section above.
CELL_CLASS_KEPT = "inh"

# Inhibitory families kept. Two groups present in the data are excluded:
# Sncg (6 reconstructed cells) and cells marked "low quality" (8 cells, no reliable family label).
# Six cells are too few to estimate a family specific baseline in the statistical models, and the
# authors themselves excluded Sncg from their classification analysis for insufficient coverage.

FAMILIES_KEPT = ["Pvalb", "Sst", "Vip", "Lamp5"]
MINIMUM_CELLS_PER_FAMILY = 20

# The repository README (errata section) states that the reconstruction for this cell is wrong, being a
# modified copy of another neuron, and recommends excluding it from follow-up analyses. It is excluded
# here. (It is an excitatory cell, so it does not in fact enter this analysis, but the filter is kept
# so that the check is explicit and reproducible.)
ERRATUM_CELLS = ["20171207_sample_1"]

# Morphology features
# Columns in the morphology table that are not shape measurements of inhibitory neurons and are
# therefore removed before the shape analysis. Columns beginning with "apical" describe the apical
# dendrite of excitatory neurons and are blank for inhibitory cells.
EXCITATORY_ONLY_PREFIX = '"apical"'
NON_SHAPE_COLUMNS = [
    "dendrite median path angle",   # blank for every inhibitory cell
    "dendrite min branch angle",    # blank for every inhibitory cell
    "normalized depth",             # position in the cortex, used as a covariate
    "normalized layer",             # position in the cortex, used as a covariate
    "cell class",                   # the inhibitory or excitatory label
]

# Number of shape axes (principal components) carried into the models.
# Chosen from the scree plot (figure panel b): the first five axes sit before the elbow and together
# explain about 55% of the shape variation; every later axis adds under 5% and would add another
# 1,500 gene tests. The variance reached by the kept axes is reported alongside.
NUMBER_OF_SHAPE_AXES = 5

# Readable names for the kept shape axes, written after inspecting the loadings reported in
# analysis_log.txt and figure panel b. An axis without an entry is named automatically from its two
# strongest loadings.
SHAPE_AXIS_NAMES = {
    "shape_axis_1": "Dendrite size and reach",
    "shape_axis_2": "Axon size and spread",
    "shape_axis_3": "Dendrite branching",
    "shape_axis_4": "Axon branching",
    "shape_axis_5": "Axon vertical extent",
}

# Gene expression
# A cell is kept only if its exon counts sum to at least this many reads. The median cell has about
# 1.3 million; two cells have under 10,000, where counts per million become dominated by sampling noise.
# The sum is taken from the count matrix itself, because the "Exon reads" column in the metadata does
# not match the matrix for every cell.
MINIMUM_EXON_READS_PER_CELL = 50_000
# A gene is kept only if it has at least one read in this many cells.
MINIMUM_CELLS_EXPRESSING_GENE = 10
# Number of highly variable genes carried into the models.
NUMBER_OF_VARIABLE_GENES = 1500
# Scale for counts per million normalisation.
COUNTS_PER_MILLION_SCALE = 1e6

# Statistics
FALSE_DISCOVERY_THRESHOLD = 0.05   # q-value below which an association is reported
NUMBER_OF_PERMUTATIONS = 200       # for the shuffled null distribution
NUMBER_OF_BOOTSTRAPS = 100         # for gene selection stability
STABILITY_THRESHOLD = 0.70         # share of bootstraps in which a gene must be selected
LABEL_CONFIDENCE_THRESHOLDS = [0.6, 0.7, 0.8]  # for the sensitivity analysis
CROSS_VALIDATION_FOLDS = 10

# Readable labels
# Raw column names in the data files mapped to the text shown on figure axes, so that every plot can be
# read without knowing the file format.
READABLE_LABELS = {
    "RNA family": "Cell family",
    "RNA type": "Transcriptomic type",
    "RNA type confidence": "Label confidence (0 to 1)",
    "normalized depth": "Depth below brain surface (0 = surface, 1 = white matter)",
    "Exon reads": "Exon reads per cell",
    "Number of genes detected": "Genes detected per cell",
    "axon total length": "Total axon length (micrometres)",
    "axon width": "Axon spread, horizontal (micrometres)",
    "axon height": "Axon spread, vertical (micrometres)",
    "axon robust width": "Axon spread, horizontal, robust (micrometres)",
    "axon robust height": "Axon spread, vertical, robust (micrometres)",
    "axon branch points": "Number of axon branch points",
    "axon tips": "Number of axon tips",
    "axon max branch order": "Maximum axon branch order",
    "axon max path distance to soma": "Longest axon path from cell body (micrometres)",
    "axon max Euclidean dist": "Furthest axon point from cell body (micrometres)",
    "axon above soma": "Share of axon above the cell body",
    "axon x-bias": "Axon left to right bias",
    "axon z-bias": "Axon up to down bias",
    "dendrite total length": "Total dendrite length (micrometres)",
    "dendrite width": "Dendrite spread, horizontal (micrometres)",
    "dendrite height": "Dendrite spread, vertical (micrometres)",
    "dendrite branch points": "Number of dendrite branch points",
    "dendrite tips": "Number of dendrite tips",
    "dendrite max path distance to soma": "Longest dendrite path from cell body (micrometres)",
    "dendrite above soma": "Share of dendrite above the cell body",
    "soma radius": "Cell body radius (micrometres)",
    "stems": "Number of primary dendrites",
    "EMD axon dendrite": "Separation of axon and dendrite (earth mover distance)",
    "Log1p fraction of axon above dendrite": "Axon above dendrite (log fraction)",
    "Log1p fraction of axon below dendrite": "Axon below dendrite (log fraction)",
}

# Colours used for the four families on every figure, so a family always looks the same from one slide to the next.
FAMILY_COLOURS = {
    "Pvalb": "#c0392b",
    "Sst": "#d68910",
    "Vip": "#7d3c98",
    "Lamp5": "#e59866",
}

def readable_label(column_name: str) -> str:
    """Return the readable label for a column, or a tidied version of its name."""
    if column_name in READABLE_LABELS:
        return READABLE_LABELS[column_name]
    return column_name.replace("_", " ").replace("log1p", "log").capitalize()