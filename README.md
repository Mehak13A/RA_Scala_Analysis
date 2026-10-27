# Which genes track the shape of a neuron?

Gene expression and morphology in inhibitory neurons of mouse motor cortex, using
the Patch-seq data of Scala et al. (2021).

Author: Mehak Agrawal.

## 1. Problem statement

Scala et al. (2021) recorded, stained and sequenced over 1,300 neurons from adult
mouse primary motor cortex, so that each neuron has a gene expression profile, a
set of electrical measurements and a traced three dimensional shape. The paper
linked gene expression to electrical behaviour in detail. Morphology received
less attention: a sparse regression from genes to shape appears only as a
supplementary figure, with no covariates.

The open question is which genes change their expression as the shape of a
neuron changes, what those genes might do, and how much any answer should be
trusted. The difficulty is a confound. Interneurons come in families (Pvalb,
Sst, Vip, Lamp5), and each family has both its own gene programme and its own
characteristic shape. A gene that correlates with shape across all cells may do
so only because it marks a family. Depth in the cortex behaves the same way.
The analysis is therefore built around one question: what do genes explain
about shape once family and position are accounted for?

## 2. Task overview

The brief asked for three things, using the released Patch-seq data:

1. Identify genes whose expression corresponds with morphological features of
   neurons, with the definition of morphology and the choice of cells left
   open.
2. Interpret the possible function of those genes and visualise the result.
3. Provide statistical support: why the chosen method suits the data, and how
   confident to be in the outputs.

Everything is delivered as one analysis script and one settings file that run
end to end from the three public data files and produce three result tables,
one figure and a log.

## 3. Repository structure

The project folder is `ra_scala_interview_2`. It contains:

`analysis.py`: the whole analysis, in numbered steps.
`config.py`: every setting and decision the analysis depends on.
`README.md`: this file.
`requirements.txt`: the Python packages and the versions used.
`data/`: the three input files from the authors' repository, listed in
Section 6.
`outputs/`: the files the analysis writes, listed in Section 8:
`results_gene_associations.csv`, `results_model_summary.csv`,
`results_function.csv`, `figure_main.png` and `analysis_log.txt`.

## 4. Scripts

`config.py` holds no computation. It contains the data and output folder paths,
the random seed, the cell selection rules (which class, which families, the
minimum cells per family, the cell flagged in the repository errata), the
morphology settings (which columns are not shape measurements, how many shape
axes to keep and their readable names), the expression settings (minimum reads
per cell, minimum cells per gene, number of variable genes), the statistical
thresholds (false discovery rate, numbers of permutations and bootstraps,
stability threshold, label confidence thresholds, cross validation folds) and
a dictionary that maps raw column names to the readable labels used on the
figure.

`analysis.py` runs seven numbered steps, each a function, and a figure
function. Step 1 selects the cells and records a funnel of how many survive
each filter. Step 2 reduces the 49 shape measurements to five named shape
axes. Step 3 normalises the gene counts and selects 1,500 highly variable
genes with the authors' own rule. Step 4 reproduces the authors' approach, a
sparse regression from genes to shape with no covariates. Step 5 is the
adjusted analysis with family, soma depth and sequencing depth as covariates,
gene by gene and for all genes at once. Step 5b tests within family
continuity at the level of individual cells. Step 6 provides the statistical
support: a permutation null, bootstrap stability and a label confidence
sensitivity analysis. Step 7 scores transcription factor activity and tests it
against shape, and attempts gene set enrichment. The main function assembles
the result tables and draws the figure. Every line printed is also written to
`analysis_log.txt`.

## 5. Dependencies

Each package and why it is needed:

pandas: tables with named rows and columns; every dataset is a DataFrame and
every result table is written with it.
numpy: fast arithmetic on arrays; used for the normalisation, the residual
calculations, the permutations and the bootstraps.
scipy: statistical functions; Pearson and Spearman correlations, the t
distribution for p values, and a smoothing spline used in the q value
estimate.
scikit-learn: the machine learning toolkit; standardisation (StandardScaler),
principal component analysis (PCA), ordinary least squares
(LinearRegression), sparse penalised regression (MultiTaskElasticNet and its
cross validated version) and the cross validation splitter (KFold).
matplotlib: draws the six panel figure.
decoupler: scores the activity of transcription factors in each cell from a
regulon database (used in step 7; needs an internet connection to download
the CollecTRI network).
gprofiler-official: queries the g:Profiler web service for gene set
enrichment (step 7; needs an internet connection).

Exact versions are in `requirements.txt`.

## 6. Data

All inputs come from the authors' public repository
https://github.com/berenslab/mini-atlas, folder `data/`. Three files are used.

`m1_patchseq_meta_data.csv` (tab separated): one row per sequenced cell with
the family and transcriptomic type assigned by the authors, the confidence of
that assignment, any exclusion reason, sequencing depth, number of genes
detected, animal, batch and Cre line.

`m1_patchseq_morph_features.csv`: one row per reconstructed cell with 63 shape
measurements extracted by the authors from the traced reconstructions, plus
normalised depth and the cell class (inhibitory or excitatory).

`m1_patchseq_exon_counts.csv.gz`: read counts per gene per cell, counting reads
that fall in exons; genes as rows and cells as columns.

Files available but not used, with the reason: the intron count file (exon
counts are the standard measure; adding introns is listed as a next step); the
depth profile file (it describes position in the cortex, which is handled as a
covariate, not shape); the electrophysiology files (the task concerns
morphology); the physiological temperature files (no reconstructed shapes); the
raw SWC reconstructions (the hosting server was unavailable during the task
window, so the authors' extracted features are used, as the panel suggested);
the raw sequencing reads (they would only reproduce the count table).

Cells used: inhibitory neurons, because they are the only class with both axon
and dendrite reconstructed; families Pvalb, Sst, Vip and Lamp5; no exclusion
flag from the authors; complete shape measurements; at least 50,000 exon reads.
The Sncg family (6 cells) and 8 unlabelled cells are excluded because six cells
cannot support a family baseline in the models; the authors excluded Sncg from
their own classifier for the same reason. The reconstruction flagged in the
repository errata is excluded (it is excitatory, so it would not have entered
anyway). The result is 358 neurons: Pvalb 144, Sst 104, Vip 63, Lamp5 47.

## 7. Methodology

### 7.1 Workflow

Step 1, cell selection. The metadata is joined onto the morphology table by
cell identifier and the filters above are applied in a fixed order, each
recorded with the number of cells remaining. Shape measurements that describe
only excitatory cells, or position rather than shape, are removed, leaving 49.
The exon counts are restricted to the kept cells and transposed to cells by
genes. Cells whose exon counts sum to under 50,000 are removed, the total being
taken from the matrix because the metadata column disagrees with it for one
cell.

Step 2, shape axes. The 49 measurements are standardised (mean zero, standard
deviation one, so that micrometres, counts and ratios contribute equally) and
reduced by principal component analysis. The first five components are kept:
they sit before the elbow of the scree plot, together explain 55 percent of
shape variation, and each has a clear anatomical reading from its loadings
(dendrite size and reach; axon size and spread; dendrite branching; axon
branching; axon vertical extent). Each axis is correlated with soma depth; the
first correlates at 0.43, which is why depth is carried as a covariate.

Step 3, expression. Counts are converted to counts per million per cell and
log transformed as log2(x + 1). Genes detected in fewer than 10 cells are
removed, leaving 20,551. The 1,500 most informative genes are selected with
the authors' rule: for each gene, the mean log count in cells where it is
clearly on and the fraction of cells where it is off are computed, and genes
above a shifted exponential curve in that plane are kept, the shift found by
binary search so that exactly 1,500 are selected. The selected genes are
standardised for the regression models.

Step 4, baseline. A multi output elastic net (a linear model predicting all
five shape axes from all 1,500 genes, with a penalty that sets most
coefficients to zero) is fitted with the penalty chosen by 10 fold cross
validation. This is the same model family as the authors' sparse reduced rank
regression. The selected genes and the cross validated R squared per axis are
recorded.

Step 5, adjusted analysis. Gene by gene: for every gene and axis a linear
model of the form axis equals gene plus family plus soma depth plus
sequencing depth is fitted and the gene coefficient, its t statistic and p
value are recorded; the 7,500 p values are converted to Storey q values. All
genes at once: three models are compared on the same cross validation folds:
family and depth only; genes only; family and depth, then genes on the
residual. The penalty for the residual model is tuned on the residuals. The
agreement between the genes only and covariates only predictions is measured
as a correlation per axis.

Step 5b, continuity. Within each family, every pair of cells is compared: the
gene distance between them (one minus the correlation of their 1,500 gene
profiles) against their shape distance (Euclidean distance on the five axes
after the covariates are removed). The correlation between the two distances
is tested with a Mantel permutation.

Step 6, support. A permutation null for the gene contribution (shape rows
shuffled within family 200 times, model refitted each time); bootstrap
stability of gene selection (100 resamples of cells); and a sensitivity
analysis repeating the gene by gene model on cells whose family label
confidence is at least 0.6, 0.7 and 0.8.

Step 7, function. Transcription factor activity per cell is scored with
decoupler on the mouse CollecTRI regulons and tested against the shape axes
with the same adjusted model and q values. Gene set enrichment with g:Profiler
is attempted on genes that are both bootstrap stable and significant.

### 7.2 Methods used, each in brief

Counts per million: each cell's gene counts divided by the cell's total and
multiplied by one million, so cells sequenced to different depths become
comparable.

Log transform, log2(x + 1): compresses the range of counts so that highly and
lowly expressed genes sit on a comparable scale; the plus one keeps zeros at
zero.

Standardisation: rescaling a variable to mean zero and standard deviation one,
so that variables measured in different units contribute equally and effect
sizes are comparable.

Principal component analysis (PCA): finds the directions in which a set of
measurements varies most, in order, each uncorrelated with the previous ones;
each cell receives a score on each direction, and the loadings say which
measurements define it.

Highly variable gene selection (the authors' rule): keeps genes that are
strongly expressed when on and often off, which are the genes that
distinguish cells; implemented as a curve in the plane of mean when on
against fraction off.

Elastic net regression (MultiTaskElasticNet): a linear regression with a
penalty combining lasso (which sets coefficients exactly to zero, giving a
sparse set of genes) and ridge (which stabilises correlated predictors); the
multi task version predicts several outputs at once and selects the same
genes for all of them.

Cross validation (10 fold): the cells are split into ten groups; the model is
fitted on nine and tested on the tenth, rotating, so that accuracy is measured
on cells the model has not seen.

R squared: the share of variation in a shape axis that a model predicts
correctly, from 0 (no better than the average) to 1 (perfect); reported on
held out cells.

Ordinary least squares with covariates: a linear model whose coefficient for
the gene is interpreted holding family, depth and sequencing depth fixed.

Frisch Waugh Lovell theorem: a result from regression theory which says that
the coefficient of one predictor in a model with several predictors can be
obtained in two stages: first remove from both the outcome and that predictor
whatever the other predictors explain (take the residuals of each after
regressing on the covariates), then regress one residual on the other. The
slope of that second regression equals the predictor's coefficient in the full
model. It is used here because it lets all 1,500 gene models for a shape axis
be fitted at once with one matrix operation on the covariates, instead of
7,500 separate regressions, while giving identical coefficients and p values.

Storey q values: a false discovery rate correction for many simultaneous
tests; q is the expected share of false findings among results at least as
strong; it also estimates pi0, the share of tests with no real effect.

Permutation test: the data are shuffled to destroy the relationship under
test, the statistic recomputed many times, and the real statistic compared to
that null distribution; the p value is the share of shuffles that match or
exceed it.

Bootstrap stability selection: the cells are resampled with replacement many
times and the sparse model refitted; a gene's stability is the share of
resamples in which it is selected.

Sensitivity analysis: repeating an analysis on a stricter subset (here cells
with high label confidence) to check that the conclusion does not depend on
uncertain inputs.

Mantel test: a correlation between two distance matrices over the same
objects (here gene distance and shape distance between pairs of cells),
tested by permuting the labels of one matrix.

Transcription factor activity scoring (decoupler, univariate linear model):
for each cell and each transcription factor, a small regression asks whether
the factor's known target genes (from the CollecTRI database) are
collectively up or down, giving an activity score per factor per cell.

Gene set enrichment (g:Profiler): tests whether a list of genes contains more
members of a known biological process than expected by chance, against a
chosen background.

## 8. Outputs

`results_gene_associations.csv`: 7,500 rows, one per gene and shape axis, with
the adjusted effect, t statistic, p value, q value, bootstrap stability,
whether the baseline model selected the gene, whether it is a known family
marker and the readable axis name. Sorted with the strongest associations
first.

`results_model_summary.csv`: five rows, one per shape axis, with the axis
description and variance explained, its correlation with depth, the cross
validated R squared of the three models, the agreement between the genes only
and covariates only predictions, the gene contribution on the residual with
its permutation null and p value, the number of genes at q below 0.05, pi0,
the label confidence sensitivity results and the within family Mantel
correlations and p values.

`results_function.csv`: the transcription factors associated with shape at q
below 0.05 (none) and the five strongest, plus any enriched gene sets (none,
because no gene list passed).

`figure_main.png`: six panels: the cells by family and depth; the five shape
axes; the R squared of the three models with permutation p; the agreement of
the two models' predictions; the 25 strongest adjusted associations as a
heatmap; the within family Mantel correlations against their nulls.

`analysis_log.txt`: everything printed during the run, including the cell
funnel, the axis descriptions, the ten strongest genes and the five strongest
transcription factors.

A detailed reading guide for every output is in `Outputs_Explanation.docx`.

## 9. How to run

1. Create a Python 3.12 environment and install the packages:
   `pip install -r requirements.txt`
2. Place the three data files from the authors' repository in the `data`
   folder.
3. Open `config.py` and set `DATA_DIR` and `OUTPUT_DIR` to the data and
   outputs folders (by default they point at `data` and `outputs` inside the
   project folder).
4. From the project folder run `python analysis.py`.

Runtime is a few minutes, most of it in the permutations. Step 7 needs an
internet connection; if it is unavailable the step is skipped with a message
and every other output is still produced. The random seed is fixed in
`config.py`, so a rerun gives the same numbers.

## 10. References

Scala F, Kobak D, Bernabucci M, Bernaerts Y, Cadwell CR, Castro JR, Hartmanis
L, Jiang X, Laturnus S, Miranda E, Mulherkar S, Tan ZH, Yao Z, Zeng H, Sandberg
R, Berens P, Tolias AS. Phenotypic variation of transcriptomic cell types in
mouse motor cortex. Nature 598, 144 to 150 (2021).
https://doi.org/10.1038/s41586-020-2907-3

Berens laboratory. mini-atlas: analysis code and preprocessed data for Scala
et al. 2021. https://github.com/berenslab/mini-atlas

Storey JD, Tibshirani R. Statistical significance for genomewide studies.
Proceedings of the National Academy of Sciences USA 100, 9440 to 9445 (2003).
