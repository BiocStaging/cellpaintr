#' cellpaintr: perturbation analysis of Cell Painting data
#'
#' cellpaintr tests whether a perturbation, such as a drug, changes the
#' morphology of cells reproducibly across a handful of replicates, such as
#' patients or plates. It learns on single cells and tests on replicates: a
#' random forest trained on the cells of all other replicates scores the
#' cells of each held-out replicate (pre-validation), and a patient-level
#' test compares the scores between the two conditions, for all features and
#' for interpretable feature groups such as measurement types or imaging
#' channels.
#'
#' @section Workflow:
#' \enumerate{
#'   \item Load a per-cell CellProfiler table with \code{\link{loadData}}.
#'   \item Remove poor cells and features with \code{\link{removeOutliers}},
#'     \code{\link{removeNAs}}, \code{\link{removeLowVariance}}, and
#'     \code{\link{removeZeroInflation}}.
#'   \item Transform and scale the features with
#'     \code{\link{transformScale}}, and optionally remove plate row and
#'     column effects with \code{\link{removeRowColEffect}}.
#'   \item Score every cell with a classifier that never saw its replicate,
#'     with \code{\link{predictLOO}}.
#'   \item Test each feature group with \code{\link{calculateStats}}: an exact
#'     permutation test that refits all models for every relabeling of the
#'     replicates, or a corrected paired t-test that needs no refits.
#'   \item Summarize the results with \code{\link{volcanoPlot}},
#'     \code{\link{plotLOO}}, \code{\link{plotROC}}, and
#'     \code{\link{plotAUC}}.
#' }
#'
#' \code{\link{generateData}} simulates a small CellProfiler table, and
#' \code{system.file("extdata", "jump_subset.csv.gz", package = "cellpaintr")}
#' holds a small subset of the JUMP Cell Painting data. The vignette walks
#' through a complete analysis:
#' \code{vignette("cellpaintr-vignette", package = "cellpaintr")}.
#'
#' @keywords internal
"_PACKAGE"
