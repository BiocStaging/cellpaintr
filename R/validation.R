#' Internal functions to validate the inputs of exported functions
#'
#' Each function returns its input invisibly or stops with an informative
#' error message from \code{cli::cli_abort()}.
#'
#' @importFrom cli cli_abort
#' @importFrom methods is
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleCellExperiment colData
#' @name validation
#' @noRd
NULL

check_sce <- function(sce, arg = "sce") {
    if (!is(sce, "SingleCellExperiment")) {
        cli_abort(c(
            "{.arg {arg}} must be a {.cls SingleCellExperiment}.",
            "x" = "You've supplied a {.cls {class(sce)}}."
        ))
    }
    invisible(sce)
}

check_string <- function(x, arg) {
    if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
        cli_abort("{.arg {arg}} must be a single non-empty string.")
    }
    invisible(x)
}

check_flag <- function(x, arg) {
    if (!is.logical(x) || length(x) != 1 || is.na(x)) {
        cli_abort("{.arg {arg}} must be {.code TRUE} or {.code FALSE}.")
    }
    invisible(x)
}

check_number <- function(x, arg, min = -Inf, integer = FALSE) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < min ||
        (integer && x != round(x))) {
        what <- if (integer) "a whole number" else "a number"
        cli_abort("{.arg {arg}} must be {what} of at least {min}.")
    }
    invisible(x)
}

check_assay <- function(sce, assay_type) {
    check_string(assay_type, "assay_type")
    if (!assay_type %in% assayNames(sce)) {
        cli_abort(c(
            "Assay {.val {assay_type}} not found.",
            "i" = "Available assays: {.val {assayNames(sce)}}."
        ))
    }
    invisible(assay_type)
}

check_column <- function(sce, column, arg) {
    check_string(column, arg)
    if (!column %in% names(colData(sce))) {
        cli_abort(c(
            "{.arg {arg}} = {.val {column}} is not a column of
            {.code colData(sce)}.",
            "i" = "Available columns: {.val {names(colData(sce))}}."
        ))
    }
    invisible(column)
}

check_levels <- function(sce, target, interest_level, reference_level) {
    check_string(interest_level, "interest_level")
    check_string(reference_level, "reference_level")
    if (identical(interest_level, reference_level)) {
        cli_abort(
            "{.arg interest_level} and {.arg reference_level} must differ."
        )
    }
    values <- as.character(sce[[target]])
    missing <- setdiff(c(interest_level, reference_level), values)
    if (length(missing) > 0) {
        cli_abort(c(
            "Not found in {.code sce${target}}: {.val {missing}}.",
            "i" = "Available levels: {.val {sort(unique(values))}}."
        ))
    }
    invisible(TRUE)
}

check_predict_args <- function(sce, target, group, interest_level,
    reference_level, types, channels, weights,
    n_threads, assay_type, progress) {
    check_sce(sce)
    check_column(sce, target, "target")
    check_column(sce, group, "group")
    if (identical(target, group)) {
        cli_abort("{.arg target} and {.arg group} must be different columns.")
    }
    check_levels(sce, target, interest_level, reference_level)
    for (arg in c("types", "channels")) {
        value <- get(arg)
        if (!is.null(value) && (!is.character(value) || anyNA(value))) {
            cli_abort("{.arg {arg}} must be {.code NULL} or a character
                vector.")
        }
    }
    if (!is.null(weights)) {
        check_column(sce, weights, "weights")
    }
    check_number(n_threads, "n_threads", min = 1, integer = TRUE)
    check_assay(sce, assay_type)
    check_flag(progress, "progress")
    invisible(TRUE)
}

check_prevalidated <- function(sce, target, group) {
    check_sce(sce)
    if (!"prevalidated" %in% SingleCellExperiment::reducedDimNames(sce)) {
        cli_abort(c(
            "{.code reducedDim(sce, \"prevalidated\")} is missing.",
            "i" = "Pass the output of {.fn predictLOO}."
        ))
    }
    check_column(sce, target, "target")
    check_column(sce, group, "group")
    invisible(TRUE)
}

check_cell_table <- function(data, cell_file, feature_vars) {
    if (!"ImageNumber" %in% names(data)) {
        cli_abort(c(
            "{.file {cell_file}} has no column {.field ImageNumber}.",
            "i" = "Pass a per-cell table exported by CellProfiler."
        ))
    }
    if (!any(startsWith(names(data), "Metadata_"))) {
        cli_abort(c(
            "{.file {cell_file}} has no {.field Metadata_} columns.",
            "i" = paste(
                "Export the image metadata with the cells (option \"Add",
                "image metadata columns to your object data file\" of",
                "ExportToSpreadsheet) or join them from the image table by",
                "{.field ImageNumber}."
            )
        ))
    }
    numeric_vars <- vapply(data[feature_vars], is.numeric, logical(1))
    if (length(feature_vars) == 0 || !all(numeric_vars)) {
        cli_abort(c(
            "Features must be numeric columns.",
            "x" = "Non-numeric: {.field {feature_vars[!numeric_vars]}}.",
            "i" = "Prefix metadata columns with {.field Metadata_}."
        ))
    }
    invisible(TRUE)
}

check_refit_settings <- function(settings, target, group) {
    if (is.null(settings) || !identical(settings$target, target) ||
        !identical(settings$group, group)) {
        cli_abort(c(
            "The permutation test refits the models of {.fn predictLOO}.",
            "i" = paste(
                "Pass its output with the same {.arg target} and",
                "{.arg group}, or use {.code method = \"t\"}."
            )
        ))
    }
    invisible(TRUE)
}

check_groups <- function(sce, group, interest_level, reference_level) {
    n_groups <- length(unique(sce[[group]]))
    if (n_groups < 2) {
        cli_abort(c(
            "{.arg group} = {.val {group}} needs at least two levels with cells
            of {.val {interest_level}} or {.val {reference_level}}.",
            "x" = "Found {n_groups}."
        ))
    }
    invisible(TRUE)
}
