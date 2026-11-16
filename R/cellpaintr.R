#' Simulate CellProfiler data and write to a temporary file
#'
#' @importFrom readr read_csv write_csv
#' @importFrom stringr str_detect
#' @importFrom stats rgamma rnorm rpois
#' @export
#'
#' @return path to csv file
#'
#' @examples
#' set.seed(23)
#' generateData()
#'
generateData <- function() {
    # use column names from CellProfiler and simulate data
    header_file <- system.file(
        "extdata", "header.csv",
        package = "cellpaintr", mustWork = TRUE
    )
    names_cp <- names(read_csv(header_file, show_col_types = FALSE))

    # create a 384-well plate
    grid <- expand.grid(col = sprintf("%02d", seq(24)), row = LETTERS[seq(16)])
    wells <- paste0(grid$row, grid$col)

    # 1,000 cells on 100 images, six patients, eight drugs, and 10 features
    # per measurement type
    n_cells <- 1000
    p <- 10
    draw <- function(x) sample(x, n_cells, replace = TRUE)
    meta <- data.frame(
        ObjectNumber = seq(n_cells),
        ImageNumber = draw(100),
        Metadata_Well = draw(wells),
        Metadata_Patient = draw(paste0("P", seq(6))),
        Metadata_Drug = draw(paste0("D", seq(8)))
    )

    # non-negative shape, intensity, and texture features
    shape <- simulate_features(names_cp, "AreaShape_", p, n_cells,
        function(n) rpois(n, lambda = 10)
    )
    intensity <- simulate_features(names_cp, "Intensity_", p, n_cells,
        function(n) rgamma(n, shape = 1)
    )
    texture <- simulate_features(names_cp, "Texture_", p, n_cells,
        function(n) rgamma(n, shape = 1)
    )
    df <- data.frame(meta, shape, intensity, texture)

    # shift the texture features of D7 and D8
    shifts <- c(D7 = 0.2, D8 = 0.4)
    for (drug in names(shifts)) {
        cells <- which(df$Metadata_Drug == drug)
        df[cells, names(texture)] <- df[cells, names(texture)] + shifts[[drug]]
    }

    # write the simulated data frame to a temporary file
    cell_file <- tempfile(fileext = ".csv")
    write_csv(df, cell_file)
    cell_file
}

#' Internal function to simulate the features of one measurement type
#'
#' @param names_cp Column names of a CellProfiler table
#' @param pattern Pattern of the feature names of the measurement type
#' @param p Number of features
#' @param n_cells Number of cells
#' @param generator Function that draws a given number of random values
#' @return \code{data.frame} with one column per sampled feature name
#' @noRd
#'
simulate_features <- function(names_cp, pattern, p, n_cells, generator) {
    features <- sample(names_cp[str_detect(names_cp, pattern)], p)
    mat <- matrix(generator(n_cells * p), nrow = n_cells, ncol = p)
    colnames(mat) <- features
    as.data.frame(mat, check.names = FALSE)
}

#' Load cell painting data from file and convert to a SingleCellExperiment
#'
#' Reads a per-cell table exported by CellProfiler's ExportToSpreadsheet
#' module, for example \code{MyExpt_Cells.csv}, with one row per cell. The
#' table needs the column \code{ImageNumber}, at least one image metadata
#' column whose name starts with \code{Metadata_} (CellProfiler adds them with
#' the option "Add image metadata columns to your object data file"), and
#' numeric measurement columns such as \code{AreaShape_Area} or
#' \code{Intensity_MeanIntensity_DNA}. Metadata columns become the column
#' data of the cells, without the \code{Metadata_} prefix, together with
#' \code{ImageNumber}, \code{ObjectNumber}, and \code{Location_Center_X},
#' \code{_Y}, and \code{_Z}; all other columns become the features in the
#' assay \code{"features"}.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom readr read_csv
#' @importFrom stringr str_detect str_remove
#' @importFrom cli cli_abort
#' @export
#'
#' @param cell_file character(1), path to a csv file from CellProfiler, e.g.,
#'   \code{MyExpt_Cells.csv}; compressed files such as \code{.csv.gz} work
#' @param ... Further arguments to \code{\link[readr]{read_csv}}, for example
#'   \code{num_threads} to limit the threads used to read large files
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'   with features as rows and cells as columns
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' loadData(cell_file)
#'
loadData <- function(cell_file, ...) {
    check_string(cell_file, "cell_file")
    if (!file.exists(cell_file)) {
        cli_abort("File {.file {cell_file}} does not exist.")
    }

    # read file with cell data and remove features if they exist
    data <- read_csv(cell_file, show_col_types = FALSE, ...)
    rmv_vars <- c("Number_Object_Number", "Parent_FilteredCell")
    data <- data[, !names(data) %in% rmv_vars]

    # split into meta data and features
    meta_cols <- str_detect(names(data), "^Metadata_")
    location <- paste0("Location_Center_", c("X", "Y", "Z"))
    meta_vars <- c(
        "ImageNumber", "ObjectNumber", names(data)[meta_cols], location
    )
    feature_vars <- names(data)[!names(data) %in% meta_vars]
    check_cell_table(data, cell_file, feature_vars)
    meta_cell <- data[, intersect(names(data), meta_vars)]
    features <- data[, feature_vars] |> t() # columns represent cells

    # add NAs to location if it doesn't exist
    for (column in setdiff(location, names(meta_cell))) {
        meta_cell[[column]] <- NA
    }

    # remove meta label from names
    names(meta_cell) <- str_remove(names(meta_cell), "^Metadata_")

    # store meta and features in bioconductor object
    SingleCellExperiment(list(features = features), colData = meta_cell)
}

#' Plot number of cells per image
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot geom_histogram scale_x_log10 ggtitle
#' @importFrom dplyr group_by tally
#' @importFrom tibble as_tibble
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param bins numeric(1), number of histogram bins
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' plotCellsPerImage(sce)
#'
plotCellsPerImage <- function(sce, bins = 100) {
    check_sce(sce)
    check_column(sce, "ImageNumber", "ImageNumber")
    check_number(bins, "bins", min = 1, integer = TRUE)
    colData(sce) |>
        as_tibble() |>
        group_by(ImageNumber) |>
        tally() |>
        ggplot(aes(n)) +
        geom_histogram(bins = bins) +
        scale_x_log10() +
        ggtitle("Number of Cells Per Image")
}

#' Plot correlations between features and principal components
#'
#' Computes the Pearson correlation of every feature with the first
#' principal components in \code{reducedDim(sce, "PCA")} and shows the
#' features most correlated with one of them.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom stats cor
#' @importFrom ggplot2 ggplot geom_tile aes scale_fill_gradient2 theme_minimal
#' @importFrom ggplot2 coord_fixed ylab
#' @importFrom dplyr mutate top_n
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param filter_by numeric(1), principal component used to select features
#' @param top numeric(1), number of features to show
#' @param pcs integer vector, principal components to show
#' @param assay_type character(1), name of the assay with the features
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce, robust = TRUE)
#' sce <- scater::runPCA(sce, exprs_values = "tfmfeatures", ncomponents = 10)
#' plotPCACor(sce, filter_by = 1)
#'
plotPCACor <- function(sce, filter_by = 1, top = 20, pcs = seq(5),
    assay_type = "tfmfeatures") {
    check_sce(sce)
    check_assay(sce, assay_type)
    if (!"PCA" %in% reducedDimNames(sce)) {
        cli_abort(c(
            "{.code reducedDim(sce, \"PCA\")} is missing.",
            "i" = "Run {.fn scater::runPCA} first."
        ))
    }
    features <- assay(sce, assay_type)
    scores <- reducedDim(sce, "PCA")[, pcs]
    features_pcs <- cor(t(as.matrix(features)), scores)

    keep <- rownames(features_pcs)
    importance <- abs(features_pcs[, filter_by])

    features_pcs <- features_pcs |>
        as_tibble() |>
        mutate(name = keep, importance = importance) |>
        top_n(top, importance)

    features_pcs <- features_pcs |>
        pivot_longer(!c(name, importance),
            names_to = "PC", values_to = "cor"
        ) |>
        mutate(PC = str_remove(PC, "PC"))

    ggplot(features_pcs, aes(PC, name, fill = cor)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(
            low = "blue", high = "red", mid = "white",
            midpoint = 0, limit = c(-1, 1), space = "Lab",
            name = "Pearson\nCorrelation"
        ) +
        theme_minimal() +
        coord_fixed() +
        ylab("feature name")
}

#' Remove cells if not enough or too many in one image
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr group_by tally filter pull
#' @importFrom tibble as_tibble
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param min numeric(1), remove the cells of images with fewer cells
#' @param max numeric(1), remove the cells of images with more cells
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- removeOutliers(sce, min = 0, max = 300)
#'
removeOutliers <- function(sce, min, max) {
    check_sce(sce)
    check_column(sce, "ImageNumber", "ImageNumber")
    check_number(min, "min")
    check_number(max, "max", min = min)
    cell_ids <-
        colData(sce) |>
        as_tibble() |>
        group_by(ImageNumber) |>
        tally() |>
        filter(n >= min & n <= max) |>
        pull(ImageNumber)
    sce[, sce$ImageNumber %in% cell_ids]
}

#' Remove cells with missing features
#'
#' @importFrom SummarizedExperiment assay
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type character(1), name of the assay with the features
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- removeNAs(sce)
#'
removeNAs <- function(sce, assay_type = "features") {
    check_sce(sce)
    check_assay(sce, assay_type)
    mat <- assay(sce, assay_type)
    sce[, colSums(is.na(mat)) == 0]
}

#' Filter low variance features
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom MatrixGenerics rowVars rowMads
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param threshold numeric(1), keep features whose spread exceeds it
#' @param robust logical(1), if \code{TRUE} the spread is the median absolute
#'   deviation, otherwise the variance
#' @param assay_type character(1), name of the assay with the features
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- removeLowVariance(sce)
#'
removeLowVariance <- function(sce, threshold = 0, robust = FALSE,
    assay_type = "features") {
    check_sce(sce)
    check_number(threshold, "threshold")
    check_flag(robust, "robust")
    check_assay(sce, assay_type)
    mat <- as.matrix(assay(sce, assay_type))
    spread <- if (robust) rowMads(mat) else rowVars(mat)
    sce[which(spread > threshold), ]
}

#' Remove zero-inflated features
#'
#' @importFrom SummarizedExperiment assay
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param proportion numeric(1) between 0 and 1, remove features with a
#'   larger proportion of zeros
#' @param assay_type character(1), name of the assay with the features
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- removeZeroInflation(sce)
#'
removeZeroInflation <- function(sce, proportion = 0.2,
                                assay_type = "features") {
    check_sce(sce)
    check_number(proportion, "proportion", min = 0)
    check_assay(sce, assay_type)
    mat <- as.matrix(assay(sce, assay_type))
    prop_zeros <- rowMeans(mat == 0, na.rm = TRUE)
    sce[!(prop_zeros > proportion), ]
}

#' Transform, center, and scale features
#'
#' Transforms each family of CellProfiler features to an approximately
#' unbounded, symmetric scale (logarithm for intensities and
#' dispersion-type texture measures, Fisher's z for correlations, logit for
#' proportions, and a rank-based inverse normal transformation for
#' thresholded colocalization coefficients), caps infinite values at the
#' finite extremes of each feature, and standardizes every feature.
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom stats qlogis qnorm
#' @importFrom MatrixGenerics rowMedians rowMads rowMaxs rowMins
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param robust logical(1), if \code{TRUE} robust z-scores (median and
#'   median absolute deviation), otherwise standard z-scores
#' @param assay_type character(1), name of the assay with the features
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'   with the transformed features in the new assay \code{"tfmfeatures"}
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
#'
transformScale <- function(sce, robust = FALSE, assay_type = "features") {
    check_sce(sce)
    check_flag(robust, "robust")
    check_assay(sce, assay_type)
    mat <- transform_families(as.matrix(assay(sce, assay_type)))
    mat <- cap_infinite(mat)
    if (!robust) {
        # option: standard z-score
        mat <- t(scale(t(mat)))
    } else {
        # option: robust z-score
        mat <- (mat - rowMedians(mat)) / rowMads(mat)
    }

    # add transformed features
    assay(sce, "tfmfeatures") <- mat
    sce
}

#' Internal function to transform each family of CellProfiler features
#'
#' @importFrom stats qlogis qnorm
#' @param mat Matrix with features as rows and cells as columns
#' @return Matrix of the same dimensions
#' @noRd
#'
transform_families <- function(mat) {
    features <- rownames(mat)

    # intensities and dispersion-type texture measures: logarithm
    log_patterns <- c(
        "Correlation_K_", "Variance_", "Contrast", "MeanIntensity",
        "StdIntensity", "IntegratedIntensity", "SumAverage",
        "MedianIntensity", "MADIntensity", "UpperQuartileIntensity",
        "LowerQuartileIntensity", "MaxIntensity", "MinIntensity"
    )
    ids <- str_detect(features, paste(log_patterns, collapse = "|"))
    mat[ids, ] <- log(pmax(mat[ids, ], 0))

    # correlation coefficients: Fisher's z-transformation
    ids <- str_detect(features, "_Correlation_")
    mat[ids, ] <- atanh(pmax(pmin(mat[ids, ], 1), -1))

    # proportions: logit
    ids <- str_detect(features, "Texture_AngularSecondMoment_")
    mat[ids, ] <- qlogis(pmax(pmin(mat[ids, ], 1), 0))

    # thresholded colocalization coefficients: rank-based inverse normal
    # transformation with Blom's offset, which avoids -Inf and Inf
    ids <- which(
        str_detect(features, "Correlation_Costes_|Correlation_Manders_")
    )
    blom <- function(x) {
        r <- rank(x, na.last = "keep", ties.method = "random")
        qnorm((r - 0.375) / (sum(!is.na(x)) + 0.25))
    }
    for (i in ids) {
        mat[i, ] <- blom(mat[i, ])
    }
    mat
}

#' Internal function to cap infinite values at the finite extremes
#'
#' @importFrom MatrixGenerics rowMaxs rowMins
#' @param mat Matrix with features as rows and cells as columns
#' @return Matrix of the same dimensions with every infinite value replaced
#'   by the largest or smallest finite value of its feature
#' @noRd
#'
cap_infinite <- function(mat) {
    finite <- mat
    finite[!is.finite(finite)] <- NA
    upper <- which(mat == Inf, arr.ind = TRUE)
    lower <- which(mat == -Inf, arr.ind = TRUE)
    mat[upper] <- rowMaxs(finite, na.rm = TRUE)[upper[, 1]]
    mat[lower] <- rowMins(finite, na.rm = TRUE)[lower[, 1]]
    mat
}

#' Remove row and column effect using the median polish procedure
#'
#' Fits Tukey's median polish to the per-well medians of each feature and
#' subtracts the overall, row, and column effects from every cell, similar
#' to the B-score of high-throughput screening. The column data need the
#' factors \code{Row} and \code{Col} with the plate position of each cell.
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom stats qlogis qnorm medpolish median
#' @importFrom scrapper aggregateAcrossCells
#' @importFrom tidyr pivot_wider
#' @importFrom cli cli_abort
#' @importFrom dplyr left_join select join_by
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param assay_type character(1), name of the assay with the features
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'   with the corrected features in the new assay \code{"corrected"}
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
#' row <- substr(sce$Well, 1, 1)
#' row <- factor(row, levels = sort(unique(row)))
#' col <- substr(sce$Well, 2, 3)
#' col <- factor(as.numeric(col))
#' sce$Row <- row
#' sce$Col <- col
#' sce <- removeRowColEffect(sce)
#'
removeRowColEffect <- function(sce,
    assay_type = "tfmfeatures") {
    # check input
    check_sce(sce)
    check_assay(sce, assay_type)
    if (!is.factor(sce$Col)) {
        cli::cli_abort(c(
            "'Col' in {.code colData(sce)$Col} must be a factor",
            "x" = "You've supplied a {.cls {class(sce$Col)}}."
        ))
    }
    if (!is.factor(sce$Row)) {
        cli::cli_abort(c(
            "'Row' in {.code colData(sce)$Row} must be a factor",
            "x" = "You've supplied a {.cls {class(sce$Row)}}."
        ))
    }

    # median of each feature per well
    mat <- assay(sce, assay_type)
    aggr <- scrapper::aggregateAcrossCells(
        mat,
        factors = colData(sce)[, c("Col", "Row")]
    )
    mat_well <- t(apply(mat, 1, function(x) tapply(x, aggr$index, median)))

    # subtract the overall, row, and column medians from every cell
    corrected <- mat
    for (i in seq_len(nrow(sce))) {
        med <- polish_wells(mat_well[i, ], aggr$combinations)
        corrected[i, ] <- mat[i, ] - med$overall -
            med$row[as.integer(sce$Row)] - med$col[as.integer(sce$Col)]
    }

    # save in sce object
    assay(sce, "corrected") <- corrected
    sce
}

#' Internal function to fit a median polish to the well medians of a feature
#'
#' @importFrom stats medpolish
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select
#' @param well_values Median of the feature in each well
#' @param combinations Row and column of each well
#' @return Output of \code{stats::medpolish()}
#' @noRd
#'
polish_wells <- function(well_values, combinations) {
    plate <- data.frame(feature = well_values, combinations) |>
        pivot_wider(names_from = Col, values_from = feature) |>
        select(-Row) |>
        as.matrix()
    # NAs occur if no cells are in a well: setting NAs to zero
    plate[is.na(plate)] <- 0
    medpolish(plate, trace.iter = FALSE)
}

#' Internal function to compute y_hat on a subset of the features
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment colData reducedDim reducedDim<-
#' @importFrom dplyr bind_cols bind_rows arrange select
#' @importFrom parsnip rand_forest set_mode set_engine
#' @importFrom rsample group_initial_split training testing
#' @importFrom purrr map
#' @importFrom ranger ranger
#' @importFrom stringr str_starts
#'
#' @param feature_name String with feature name
#' @param sce_feature \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @param starts Starting string
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param group character(1), name of the column of \code{colData(sce)}
#'   with the replicate of each cell, such as the patient; the replicates
#'   are the folds of the leave-one-out cross-validation
#' @param interest_level character(1), level of `target` compared with the
#'   reference, such as a drug
#' @param reference_level character(1), level of `target` used as the
#'   reference, such as the control
#' @param weights \code{NULL} or character(1), name of the column of
#'   \code{colData(sce)} with case weights, for example cell counts when
#'   cells are aggregated to images
#' @param n_threads numeric(1), number of threads of each random forest fit
#' @param assay_type character(1), name of the assay with the features
#' @return \code{\link[tibble]{tibble}} data frame
#' @noRd
#'
compute_y_hat <- function(feature_name,
    sce_feature,
    starts,
    target,
    group,
    interest_level,
    reference_level,
    weights,
    n_threads,
    assay_type) {
    ml_data <- group_data(feature_name, sce_feature, starts, target,
        interest_level, reference_level, assay_type
    )

    # compute y hats
    patients <- unique(sce_feature[[group]])
    result_list <- lapply(patients, function(patient) {
        train_ids <- which(sce_feature[[group]] != patient)
        test_ids <- which(sce_feature[[group]] == patient)
        if (!is.null(weights)) {
            case_weights <- sce_feature[[weights]][train_ids]
        } else {
            case_weights <- NULL
        }
        data.frame(
            cell = test_ids,
            score = fit_predict(
                ml_data, train_ids, test_ids, case_weights, n_threads
            )
        )
    })

    # combine folds
    result <- result_list |>
        bind_rows() |>
        arrange(cell)

    # rename to feature type / channel
    names(result)[2] <- feature_name
    result |> select(all_of(feature_name))
}

#' Internal function to prepare the data of one feature group
#'
#' @param feature_name String with feature name, or "all" for all features
#' @param sce_feature \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @param starts Select features that start with `feature_name` (types) or
#'   contain it (channels)
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param interest_level character(1), level of `target` compared with the
#'   reference, such as a drug
#' @param reference_level character(1), level of `target` used as the
#'   reference, such as the control
#' @param assay_type character(1), name of the assay with the features
#' @return \code{data.frame} with the factor `Target` and one column per feature
#' @noRd
#'
group_data <- function(feature_name,
    sce_feature,
    starts,
    target,
    interest_level,
    reference_level,
    assay_type) {
    # convert to regular expression
    if (feature_name == "all") {
        pattern <- ".*"
    } else {
        pattern <- feature_name
    }

    # subset features
    features <- names(sce_feature)
    if (starts) {
        features <- features[str_starts(features, pattern)]
    } else {
        features <- features[str_detect(features, pattern)]
    }
    sce_feature <- sce_feature[features, ]

    # prepare data frame with target variable and predictor matrix
    X <- assay(sce_feature, assay_type) |> t()
    cells <- colData(sce_feature) |> as_tibble()
    ml_data <- data.frame(Target = cells |> pull(target), X)
    ml_data$Target <- factor(ml_data$Target,
        levels = c(interest_level, reference_level)
    )
    ml_data
}

#' Internal function to fit a random forest and score held-out cells
#'
#' @param ml_data Output of \code{group_data}
#' @param train_ids Rows of `ml_data` to train on
#' @param test_ids Rows of `ml_data` to score
#' @param case_weights Weights of the training rows, or `NULL`
#' @param n_threads numeric(1), number of threads of each random forest fit
#' @importFrom stats predict
#' @return Numeric vector with the predicted probability of the interest level
#'   for each row in `test_ids`
#' @noRd
#'
fit_predict <- function(ml_data, train_ids, test_ids, case_weights, n_threads) {
    rf_spec <-
        parsnip::rand_forest() |>
        parsnip::set_mode("classification") |>
        parsnip::set_engine("ranger",
            num.threads = n_threads,
            case.weights = case_weights
        )
    rf_fit <- rf_spec |>
        parsnip::fit(Target ~ ., data = ml_data[train_ids, ])

    # the first column is the interest level, the first factor level
    predict(rf_fit, ml_data[test_ids, ], type = "prob")[[1]]
}

#' Predict target from features
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment colData reducedDim reducedDim<-
#' @importFrom S4Vectors metadata<-
#' @importFrom dplyr bind_cols bind_rows arrange select
#' @importFrom parsnip rand_forest set_mode set_engine
#' @importFrom rsample group_initial_split training testing
#' @importFrom purrr map
#' @importFrom ranger ranger
#' @importFrom stringr str_starts
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param group character(1), name of the column of \code{colData(sce)}
#'   with the replicate of each cell, such as the patient; the replicates
#'   are the folds of the leave-one-out cross-validation
#' @param interest_level character(1), level of `target` compared with the
#'   reference, such as a drug
#' @param reference_level character(1), level of `target` used as the
#'   reference, such as the control
#' @param types \code{NULL} or character vector of CellProfiler measurement
#'   types, such as \code{"Intensity"}; each type is a feature group of the
#'   features whose names start with it
#' @param channels \code{NULL} or character vector of imaging channels, such
#'   as \code{"DNA"} or \code{"Mito"}; each channel is a feature group of the
#'   features whose names contain it
#' @param weights \code{NULL} or character(1), name of the column of
#'   \code{colData(sce)} with case weights, for example cell counts when
#'   cells are aggregated to images
#' @param n_threads numeric(1), number of threads of each random forest fit
#' @param assay_type character(1), name of the assay with the features
#' @param progress logical(1), show a progress bar while fitting the models
#'   of feature types and channels
#' @return \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'   with the cells of the two levels, the leave-one-out prediction scores in
#'   \code{reducedDim(sce, "prevalidated")}, and the settings of this call in
#'   \code{metadata(sce)$predictLOO}, which \code{\link{calculateStats}} uses
#'   to refit the models in its permutation test
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
predictLOO <- function(sce,
    target,
    group,
    interest_level,
    reference_level,
    types = NULL,
    channels = NULL,
    weights = NULL,
    n_threads = 1,
    assay_type = "tfmfeatures",
    progress = TRUE) {
    check_predict_args(
        sce, target, group, interest_level, reference_level, types,
        channels, weights, n_threads, assay_type, progress
    )

    # subset for binary classification
    sce_subset <- sce[, sce[[target]] %in% c(reference_level, interest_level)]

    # remove unused levels; a character target becomes a factor
    sce_subset[[target]] <- droplevels(as.factor(sce_subset[[target]]))
    check_groups(sce_subset, group, interest_level, reference_level)

    # leave-one-out scores of all features, the types, and the channels;
    # types select features by prefix, channels by any match
    groups <- c("all", types, channels)
    starts <- c(TRUE, rep(TRUE, length(types)), rep(FALSE, length(channels)))
    y_hat <- purrr::map2(groups, starts, function(feature_name, starts) {
        compute_y_hat(feature_name, sce_subset, starts, target, group,
            interest_level, reference_level, weights, n_threads, assay_type
        )
    }, .progress = progress && length(groups) > 1) |>
        bind_cols() |>
        as.matrix()

    # add y_hat to sce object, keeping the cell names if there are any
    rownames(y_hat) <- colnames(sce_subset)
    reducedDim(sce_subset, "prevalidated") <- y_hat

    # record the settings, which calculateStats() needs to refit the models
    metadata(sce_subset)$predictLOO <- mget(c(
        "target", "group", "interest_level", "reference_level", "types",
        "channels", "weights", "n_threads", "assay_type"
    ))
    sce_subset
}

#' Aggregate predicted leave-one-out probabilities over meta variables
#'
#' @importFrom scrapper aggregateAcrossCells
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param group character(1), name of the column of \code{colData(sce)}
#'   with the replicate of each cell, such as the patient; the replicates
#'   are the folds of the leave-one-out cross-validation
#' @param assay_type character(1), name of the assay with the features
#' @return \code{data.frame}
#' @noRd
#'
aggregateYhat <- function(sce,
    target,
    group,
    assay_type = "tfmfeatures") {
    meta_vars <- c(target, group)

    # get group index
    aggr <- scrapper::aggregateAcrossCells(
        assay(sce, assay_type),
        factors = colData(sce)[, meta_vars]
    )

    # calculate mean score
    idx <- aggr$index
    pre <- reducedDim(sce, "prevalidated")
    agg_pre <- rowsum(pre, idx) / tabulate(idx)

    # combine it with meta vars
    cbind(
        aggr$combinations,
        agg_pre
    ) |> data.frame()
}

#' Plot predicted leave-one-out probabilities
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter ylab facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param group character(1), name of the column of \code{colData(sce)}
#'   with the replicate of each cell, such as the patient; the replicates
#'   are the folds of the leave-one-out cross-validation
#' @param assay_type character(1), name of the assay with the features
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
#' plotLOO(sce_single, target = "Drug", group = "Patient")
#'
plotLOO <- function(sce,
                    target,
                    group,
                    assay_type = "tfmfeatures") {
    check_prevalidated(sce, target, group)
    check_assay(sce, assay_type)
    meta_vars <- c(target, group)
    y_hat <- aggregateYhat(sce, target, group, assay_type)

    y_hat |>
        pivot_longer(
            cols = -c(all_of(meta_vars)),
            names_to = "features", values_to = "value"
        ) |>
        ggplot(aes(.data[[target]], value, color = .data[[target]])) +
        geom_boxplot(width = 0.2, outliers = FALSE) +
        geom_jitter(width = 0.1) +
        ylab("predicted leave-one-out probability") +
        facet_wrap(~features)
}

#' Internal function to compute the test statistic of each feature group
#'
#' One-sided t-statistic for higher mean scores in the interest level,
#' computed on the scores averaged per patient and level: paired if every
#' patient has both levels, with pooled variance otherwise.
#'
#' @importFrom dplyr pull bind_rows all_of select
#' @importFrom tidyr pivot_wider
#' @importFrom stats t.test
#'
#' @param sce Output of \code{\link{predictLOO}}
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param group character(1), name of the column of \code{colData(sce)}
#'   with the replicate of each cell, such as the patient; the replicates
#'   are the folds of the leave-one-out cross-validation
#' @param assay_type character(1), name of the assay with the features
#' @param interest_level character(1), level of `target` compared with the
#'   reference, such as a drug
#' @param reference_level character(1), level of `target` used as the
#'   reference, such as the control
#' @return \code{data.frame} with one row per feature group
#' @noRd
#'
score_statistics <- function(sce,
    target,
    group,
    assay_type,
    interest_level,
    reference_level) {
    meta_vars <- c(target, group)
    y_hat <- aggregateYhat(sce, target, group, assay_type)
    feature_vars <- setdiff(names(y_hat), meta_vars)

    lapply(feature_vars, function(feature) {
        wide <- y_hat |>
            select(all_of(c(meta_vars, feature))) |>
            pivot_wider(
                names_from = all_of(target),
                values_from = all_of(feature)
            )

        x <- pull(wide, all_of(reference_level))
        y <- pull(wide, all_of(interest_level))
        paired <- sum(is.na(c(x, y))) == 0
        result <- tryCatch(
            t.test(x, y,
                paired = paired, var.equal = TRUE,
                alternative = "less"
            ),
            error = function(e) list(statistic = NA_real_, p.value = NA_real_)
        )

        data.frame(
            Feature = feature,
            statistic = -unname(result$statistic),
            pvalue = result$p.value,
            log2FoldChange = log2(mean(y, na.rm = TRUE) / mean(x, na.rm = TRUE))
        )
    }) |> bind_rows()
}

#' Internal function to list relabelings of whole patients
#'
#' If `symmetric`, swapping the levels of all patients gives another valid
#' relabeling. It turns every score into one minus the score and exchanges
#' the two levels, so it leaves the statistic unchanged for a learner that
#' treats the two levels alike. The relabelings then come in pairs with the
#' same statistic, and only those that keep the labels of one patient are
#' listed: the first patient in a paired or unpaired design, and the first
#' patient with both levels in a mixed design. They form a subgroup of all
#' relabelings, so the test stays exact for any learner, at half the cost.
#'
#' @importFrom utils combn
#'
#' @param design Either "paired", "unpaired", or "mixed"
#' @param n Number of patients
#' @param observed Logical vector of patients with the interest level
#' @param exact Enumerate all relabelings instead of sampling them
#' @param n_permutations Number of random relabelings if not exact
#' @param both Logical vector of patients with both levels (mixed design)
#' @param symmetric Keep the labels of one patient fixed
#' @return Logical matrix with one row per relabeling, the observed labels
#'   excluded, and one column per patient: patients whose two levels are
#'   swapped (paired design, and patients with both levels in a mixed
#'   design) or patients with the interest level (unpaired design, and
#'   patients with one level in a mixed design)
#' @noRd
#'
patient_relabelings <- function(design, n, observed, exact, n_permutations,
                                both = NULL, symmetric = FALSE) {
    if (design == "mixed") {
        return(mixed_relabelings(
            observed, both, exact, n_permutations, symmetric
        ))
    }
    if (design == "paired") {
        # swaps of the free patients; the first one keeps its labels if
        # symmetric
        free <- if (symmetric) seq_len(n)[-1] else seq_len(n)
        grid <- matrix(FALSE,
            nrow = if (exact) 2^length(free) else n_permutations, ncol = n
        )
        grid[, free] <- if (exact) {
            as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(free))))
        } else {
            sample(c(FALSE, TRUE), n_permutations * length(free),
                replace = TRUE
            )
        }
        if (exact) {
            grid <- grid[-1, , drop = FALSE]
        }
        return(unname(grid))
    }
    # unpaired: the first patient keeps its level if symmetric
    free <- if (symmetric) seq_len(n)[-1] else seq_len(n)
    n_free_interest <- sum(observed[free])
    grid <- matrix(observed, nrow = 1, ncol = n)
    if (exact) {
        assignments <- all_subsets(length(free), n_free_interest)
        grid <- grid[rep(1, nrow(assignments)), , drop = FALSE]
        grid[, free] <- assignments
        is_observed <- apply(grid, 1, function(row) all(row == observed))
        return(grid[!is_observed, , drop = FALSE])
    }
    grid <- grid[rep(1, n_permutations), , drop = FALSE]
    grid[, free] <- matrix(
        replicate(
            n_permutations,
            seq_along(free) %in% sample(length(free), n_free_interest)
        ),
        ncol = length(free), byrow = TRUE
    )
    grid
}

#' Internal function to list all subsets of k out of m patients
#'
#' @importFrom utils combn
#' @return Logical matrix with one row per subset and one column per patient
#' @noRd
#'
all_subsets <- function(m, k) {
    if (k == 0 || k == m) {
        return(matrix(rep(k == m, m), nrow = 1))
    }
    t(apply(combn(m, k), 2, function(subset) seq_len(m) %in% subset))
}

#' Internal function to list relabelings in a mixed design
#'
#' Patients with both levels keep or swap their two levels, and the levels
#' of the patients with one level are permuted among them.
#'
#' @param observed Logical vector of patients with the interest level
#' @param both Logical vector of patients with both levels
#' @param exact Enumerate all relabelings instead of sampling them
#' @param n_permutations Number of random relabelings if not exact
#' @param symmetric The first patient with both levels keeps its labels
#' @return Logical matrix as returned by \code{patient_relabelings}
#' @noRd
#'
mixed_relabelings <- function(observed, both, exact, n_permutations,
    symmetric = FALSE) {
    single <- !both
    # patients with both levels that may swap them
    free <- both
    if (symmetric) {
        free[which(both)[1]] <- FALSE
    }
    n_free <- sum(free)
    n_single <- sum(single)
    n_interest <- sum(observed[single])
    if (!exact) {
        grid <- matrix(FALSE, nrow = n_permutations, ncol = length(both))
        grid[, free] <- sample(c(FALSE, TRUE), n_permutations * n_free,
            replace = TRUE
        )
        grid[, single] <- matrix(
            replicate(
                n_permutations,
                seq_len(n_single) %in% sample(n_single, n_interest)
            ),
            ncol = n_single, byrow = TRUE
        )
        return(grid)
    }
    swaps <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n_free)))
    assignments <- all_subsets(n_single, n_interest)
    index <- expand.grid(
        swap = seq_len(nrow(swaps)),
        assignment = seq_len(nrow(assignments))
    )
    grid <- matrix(FALSE, nrow = nrow(index), ncol = length(both))
    grid[, free] <- swaps[index$swap, , drop = FALSE]
    grid[, single] <- assignments[index$assignment, , drop = FALSE]
    is_observed <- !apply(grid[, both, drop = FALSE], 1, any) &
        apply(grid[, single, drop = FALSE], 1, function(row) {
            all(row == observed[single])
        })
    grid[!is_observed, , drop = FALSE]
}

#' Internal function to describe the design of the patients
#'
#' Paired: every patient has both levels; unpaired: every patient has one;
#' mixed: some patients have both levels, the others have one. If swapping
#' the levels of all patients is a valid relabeling, it leaves the statistic
#' unchanged, so one patient keeps its labels (`symmetric`), which halves the
#' relabelings without changing the test.
#'
#' @param sce Output of \code{\link{predictLOO}}
#' @param settings Settings of \code{\link{predictLOO}} from its metadata
#' @return List with the labels and patients of the cells, the sorted patient
#'   ids, the design, the total number of relabelings, and the symmetry
#' @noRd
#'
relabeling_design <- function(sce, settings) {
    labels <- as.character(sce[[settings$target]])
    patients <- as.character(sce[[settings$group]])
    ids <- sort(unique(patients))
    n <- length(ids)
    in_interest <- ids %in% patients[labels == settings$interest_level]
    in_reference <- ids %in% patients[labels == settings$reference_level]
    both <- in_interest & in_reference
    if (all(both)) {
        design <- "paired"
        n_total <- 2^n
    } else if (all(xor(in_interest, in_reference))) {
        design <- "unpaired"
        n_total <- choose(n, sum(in_interest))
    } else {
        design <- "mixed"
        n_total <- 2^sum(both) * choose(sum(!both), sum(in_interest & !both))
    }
    single <- !both
    symmetric <- 2 * sum(in_interest & single) == sum(single)
    if (symmetric) {
        n_total <- n_total / 2
    }
    list(
        labels = labels, patients = patients, ids = ids, n = n,
        in_interest = in_interest, both = both, design = design,
        n_total = n_total, symmetric = symmetric,
        fixed = if (design == "mixed") which(both)[1] else 1
    )
}

#' Internal function to report the relabelings of the permutation test
#'
#' @param d Output of \code{relabeling_design}
#' @param exact Whether all relabelings are enumerated
#' @param n_relabelings Number of relabelings, the observed one included
#' @param settings Settings of \code{\link{predictLOO}} from its metadata
#' @return `NULL`, invisibly; called for its messages
#' @noRd
#'
report_relabelings <- function(d, exact, n_relabelings, settings) {
    group <- settings$group
    target <- settings$target
    message(
        "Permutation test with ", if (exact) "all " else "", n_relabelings,
        " relabelings of ", d$n, " levels of '", group, "' (", d$design,
        " design); smallest attainable p-value 1/", n_relabelings, "."
    )
    if (d$symmetric) {
        message(
            "'", d$ids[d$fixed], "' keeps its labels, because exchanging the ",
            "two levels of '", target, "' in all levels of '", group,
            "' at once gives the same statistic."
        )
    }
    if (d$design == "mixed") {
        message(
            "Levels of '", group, "' with cells of only one level of '",
            target, "': ", paste(d$ids[!d$both], collapse = ", "), "."
        )
    }
    invisible(NULL)
}

#' Internal function to relabel the cells of whole patients
#'
#' @param row Logical vector, one row of \code{patient_relabelings}
#' @param d Output of \code{relabeling_design}
#' @param interest_level character(1), level of `target` compared with the
#'   reference, such as a drug
#' @param reference_level character(1), level of `target` used as the
#'   reference, such as the control
#' @return Character vector with the new label of each cell
#' @noRd
#'
relabel_cells <- function(row, d, interest_level, reference_level) {
    patients <- d$patients
    ids <- d$ids
    if (d$design == "unpaired") {
        return(ifelse(patients %in% ids[row], interest_level, reference_level))
    }
    # swap the two levels of the selected patients with both levels
    swap <- if (d$design == "paired") row else row & d$both
    swapped <- d$labels
    selected <- patients %in% ids[swap]
    swapped[selected] <- ifelse(d$labels[selected] == interest_level,
        reference_level, interest_level
    )
    if (d$design == "mixed") {
        # permute the levels among the patients with one level
        single <- patients %in% ids[!d$both]
        swapped[single] <- ifelse(patients[single] %in% ids[row & !d$both],
            interest_level, reference_level
        )
    }
    swapped
}

#' Internal function to refit all models for one relabeling
#'
#' @param row Logical vector, one row of \code{patient_relabelings}
#' @param sce Output of \code{\link{predictLOO}}
#' @param settings Settings of \code{\link{predictLOO}} from its metadata
#' @param d Output of \code{relabeling_design}
#' @param features Feature groups in the order of the observed statistics
#' @return Numeric vector with the statistic of each feature group
#' @noRd
#'
refit_statistics <- function(row, sce, settings, d, features) {
    relabeled <- sce
    relabeled[[settings$target]] <- factor(
        relabel_cells(
            row, d, settings$interest_level, settings$reference_level
        ),
        levels = levels(sce[[settings$target]])
    )
    refit <- do.call(
        predictLOO,
        c(list(sce = relabeled), settings, list(progress = FALSE))
    )
    statistics <- score_statistics(refit, settings$target, settings$group,
        settings$assay_type, settings$interest_level, settings$reference_level
    )
    statistics$statistic[match(features, statistics$Feature)]
}

#' Internal function to compute permutation p-values by relabeling patients
#'
#' The relabelings are refitted with \code{future.apply::future_lapply()},
#' sequentially unless the user sets a parallel plan with
#' \code{future::plan()}.
#'
#' @importFrom future.apply future_lapply
#'
#' @param sce Output of \code{\link{predictLOO}}
#' @param observed Output of \code{score_statistics} for the observed labels
#' @param settings Settings of \code{\link{predictLOO}} from its metadata
#' @param n_permutations Maximal number of relabelings besides the observed
#' @return List with p-values and the number of relabelings
#' @noRd
#'
permutation_pvalues <- function(sce, observed, settings, n_permutations) {
    d <- relabeling_design(sce, settings)
    exact <- d$n_total <= n_permutations + 1
    relabelings <- patient_relabelings(
        d$design, d$n, d$in_interest, exact, n_permutations, d$both,
        d$symmetric
    )
    n_relabelings <- nrow(relabelings) + 1
    report_relabelings(d, exact, n_relabelings, settings)

    # refit all leave-one-out models for each relabeling
    null_statistics <- future_lapply(
        seq_len(nrow(relabelings)),
        function(b) {
            refit_statistics(relabelings[b, ], sce, settings, d,
                observed$Feature
            )
        },
        future.seed = TRUE
    )
    null_statistics <- matrix(unlist(null_statistics),
        nrow = nrow(relabelings), byrow = TRUE
    )

    # a relabeling counts against the observed labels if its statistic is at
    # least as large; undefined statistics count against them too
    exceeds <- sweep(null_statistics, 2, observed$statistic, ">=") |
        is.na(null_statistics)
    list(
        pvalue = (1 + colSums(exceeds)) / n_relabelings,
        relabelings = n_relabelings
    )
}

#' Internal function for the factor by which the t-statistic is divided
#'
#' `"corrected"`: if each score difference is a sum of pairwise interactions
#' with the other patients, \eqn{d_p = \sum_{q \neq p} a_{pq}} with symmetric
#' \eqn{a_{pq} = a_{qp}} uncorrelated across pairs and of mean zero, then two
#' score differences have correlation \eqn{1/(n-1)}, and the t-statistic is
#' too large by the factor \eqn{\sqrt{2(n-1)/(n-2)}}.
#'
#' `"nadeau_bengio"`: the corrected resampled t-test of Nadeau and Bengio
#' (2003) multiplies the variance by \eqn{1 + J n_2 / n_1}; with \eqn{J = n}
#' splits and one held-out patient per split, the t-statistic is divided by
#' \eqn{\sqrt{(2n-1)/(n-1)}}, which assumes a correlation of \eqn{1/n}.
#'
#' @param n Number of patients
#' @param method Either `"corrected"` (at least three patients) or
#'   `"nadeau_bengio"` (at least two patients)
#' @return The factor by which the t-statistic is divided
#' @noRd
#'
correction_factor <- function(n, method = c("corrected", "nadeau_bengio")) {
    method <- match.arg(method)
    switch(method,
        corrected = sqrt(2 * (n - 1) / (n - 2)),
        nadeau_bengio = sqrt((2 * n - 1) / (n - 1))
    )
}

#' Internal function to compute p-values of a t-test with a correction factor
#'
#' @importFrom cli cli_abort
#' @importFrom stats pt
#'
#' @param sce Output of \code{\link{predictLOO}}
#' @param observed Output of \code{score_statistics}
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param group character(1), name of the column of \code{colData(sce)}
#'   with the replicate of each cell, such as the patient; the replicates
#'   are the folds of the leave-one-out cross-validation
#' @param interest_level character(1), level of `target` compared with the
#'   reference, such as a drug
#' @param reference_level character(1), level of `target` used as the
#'   reference, such as the control
#' @param method Either `"corrected"` or `"nadeau_bengio"`
#' @return Vector of p-values, one per row of `observed`
#' @noRd
#'
corrected_pvalues <- function(sce, observed, target, group,
    interest_level, reference_level,
    method = "corrected") {
    name <- if (method == "corrected") {
        "The corrected t-test"
    } else {
        "The Nadeau-Bengio t-test"
    }
    min_n <- if (method == "corrected") 3 else 2
    labels <- as.character(sce[[target]])
    patients <- as.character(sce[[group]])
    in_levels <- labels %in% c(interest_level, reference_level)
    ids <- sort(unique(patients[in_levels]))
    n <- length(ids)

    paired <- vapply(ids, function(id) {
        all(c(interest_level, reference_level) %in% labels[patients == id])
    }, logical(1))
    if (!all(paired)) {
        cli_abort(c(
            paste(
                "{name} needs a paired design, in which every",
                "level of {.arg group} has cells of both levels of",
                "{.arg target}."
            ),
            "i" = "Levels with cells of one level only: {ids[!paired]}.",
            "i" = "Use {.code method = \"permutation\"} instead."
        ))
    }
    if (n < min_n) {
        cli_abort("{name} needs at least {min_n} levels of {.arg group}.")
    }
    pt(observed$statistic / correction_factor(n, method), n - 1,
        lower.tail = FALSE
    )
}

#' Internal function for the two levels compared by calculateStats
#'
#' The levels passed to \code{\link{predictLOO}}; for objects without its
#' settings, the first two factor levels, reference first.
#'
#' @param sce Output of \code{\link{predictLOO}}
#' @param settings Settings of \code{\link{predictLOO}} from its metadata
#' @param target Name of the target column
#' @return List with `interest` and `reference`
#' @noRd
#'
analysis_levels <- function(sce, settings, target) {
    if (!is.null(settings) && identical(settings$target, target)) {
        return(list(
            interest = settings$interest_level,
            reference = settings$reference_level
        ))
    }
    list(
        interest = levels(sce[[target]])[2],
        reference = levels(sce[[target]])[1]
    )
}

#' Test each feature group for a difference between two levels
#'
#' Average the leave-one-out prediction scores of \code{\link{predictLOO}}
#' per patient and level, and test in each feature group whether the
#' interest level has higher scores than the reference level.
#'
#' The test statistic is the one-sided t-statistic for higher mean scores in
#' the interest level. It is paired if every patient (`group`) has cells of
#' both levels, and uses a pooled variance otherwise.
#'
#' With `method = "permutation"` (the default), the p-value comes from
#' relabeling whole patients: in a paired design, the two levels of a patient
#' are swapped; in an unpaired design, the levels are permuted across
#' patients; and in a mixed design, where only some patients have cells of
#' both levels, the levels of those patients are swapped and the levels of
#' the other patients are permuted among them. For each relabeling, all
#' leave-one-out models are refitted with
#' the settings of the call to \code{\link{predictLOO}}, and the p-value is
#' the proportion of relabelings, the observed one included, with a statistic
#' at least as large. Exchanging the two levels in all patients at once turns
#' every score into one minus the score and leaves the statistic unchanged,
#' so relabelings come in pairs with the same statistic. Whenever this
#' exchange is itself a relabeling (always in a paired design, and if the
#' patients with one level are split evenly between the two levels
#' otherwise), one patient keeps its labels, which halves the number of
#' relabelings and refits without changing the test. If there are at most
#' `n_permutations + 1` relabelings (\eqn{2^{n-1}} for \eqn{n} patients in a
#' paired design), all of them are enumerated and the test is exact;
#' otherwise, `n_permutations` relabelings are drawn at random. The smallest
#' attainable p-value is one over the number of relabelings, \eqn{2^{1-n}}
#' in a paired design. Because all models are refitted for every
#' relabeling, the test takes that many times longer than
#' \code{\link{predictLOO}}. Use \code{\link[base]{set.seed}} for
#' reproducible results.
#'
#' With `method = "corrected"`, the p-value comes from a corrected paired
#' t-test, which needs a paired design with at least three patients. Each
#' patient is in the training set of every other patient's model, which makes
#' the score differences of different patients positively correlated, so the
#' t-statistic is too large. If each score difference is a sum of pairwise
#' interactions with the other patients, two score differences have
#' correlation \eqn{1/(n-1)} for \eqn{n} patients, and the t-statistic is too
#' large by the factor \eqn{\sqrt{2(n-1)/(n-2)}}. The corrected t-test
#' divides the t-statistic by this factor and refers it to the t-distribution
#' with \eqn{n-1} degrees of freedom. It needs no refits, and its p-values
#' are not limited to multiples of \eqn{2^{1-n}}, but it is approximate.
#'
#' With `method = "nadeau_bengio"`, the t-statistic is divided instead by
#' \eqn{\sqrt{(2n-1)/(n-1)}}, the factor of the corrected resampled t-test of
#' Nadeau and Bengio (2003, Machine Learning 52:239-281) for \eqn{n} splits
#' with one held-out patient each, which assumes a correlation of
#' \eqn{1/n} instead of \eqn{1/(n-1)}. It corrects less and therefore has
#' more power, which pays off if part of each score difference is independent
#' of the other patients, for example noise from the random seeds of the
#' forests. If the score differences consist mostly of pairwise interactions,
#' it can exceed the nominal level with four to six patients, where
#' `method = "corrected"` holds it. It needs a paired design with at least
#' two patients and no refits.
#'
#' With `method = "t"`, the p-value comes from the t-distribution. This is
#' fast but anti-conservative, because the score differences of different
#' patients are positively correlated.
#'
#' @importFrom S4Vectors metadata
#' @importFrom cli cli_abort
#' @export
#'
#' @param sce Output of \code{\link{predictLOO}}
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param group character(1), name of the column of \code{colData(sce)}
#'   with the replicate of each cell, such as the patient; the replicates
#'   are the folds of the leave-one-out cross-validation
#' @param assay_type character(1), name of the assay with the features
#' @param method Either `"permutation"` (default) for a permutation test that
#'   relabels whole patients, `"corrected"` for the corrected paired t-test,
#'   `"nadeau_bengio"` for the paired t-test with the correction of Nadeau
#'   and Bengio, or `"t"` for p-values from the t-distribution
#' @param n_permutations Largest number of relabelings besides the observed
#'   labels; if there are more, this many are drawn at random
#' @return \code{data.frame} with one row per feature group and the columns
#'   `Target` (interest level), `Feature` (feature group), `pvalue`, and
#'   `log2FoldChange` (log2 ratio of the mean scores of the interest and the
#'   reference level). With the permutation test, the attribute
#'   `relabelings` holds the number of relabelings, the observed one
#'   included.
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     n_threads = 1
#' )
#'
#' # permutation test; by default, it enumerates all 64 relabelings of the
#' # six patients, here 9 random relabelings keep the example fast
#' calculateStats(sce_single,
#'     target = "Drug", group = "Patient",
#'     n_permutations = 9
#' )
#'
#' # corrected t-test: no refits
#' calculateStats(sce_single,
#'     target = "Drug", group = "Patient",
#'     method = "corrected"
#' )
#'
#' # correction of Nadeau and Bengio: no refits, corrects less
#' calculateStats(sce_single,
#'     target = "Drug", group = "Patient",
#'     method = "nadeau_bengio"
#' )
#'
#' # p-values from the t-distribution: fast but anti-conservative
#' calculateStats(sce_single,
#'     target = "Drug", group = "Patient",
#'     method = "t"
#' )
#'
calculateStats <- function(sce,
    target,
    group,
    assay_type = "tfmfeatures",
    method = c("permutation", "corrected", "nadeau_bengio", "t"),
    n_permutations = 999) {
    method <- match.arg(method)
    check_prevalidated(sce, target, group)
    check_assay(sce, assay_type)
    settings <- metadata(sce)$predictLOO
    levels <- analysis_levels(sce, settings, target)
    interest_level <- levels$interest
    reference_level <- levels$reference

    observed <- score_statistics(sce, target, group, assay_type,
        interest_level, reference_level
    )
    stats <- data.frame(
        Target = interest_level,
        Feature = observed$Feature,
        pvalue = observed$pvalue,
        log2FoldChange = observed$log2FoldChange
    )
    if (method == "t") {
        return(stats)
    }
    if (method %in% c("corrected", "nadeau_bengio")) {
        stats$pvalue <- corrected_pvalues(sce, observed, target, group,
            interest_level, reference_level, method
        )
        return(stats)
    }

    check_refit_settings(settings, target, group)
    check_number(n_permutations, "n_permutations", min = 1, integer = TRUE)
    permutation <- permutation_pvalues(sce, observed, settings, n_permutations)
    stats$pvalue <- permutation$pvalue
    attr(stats, "relabelings") <- permutation$relabelings
    stats
}

#' Plot predicted leave-one-out probabilities
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_hline xlab ylab geom_point
#' @importFrom dplyr mutate
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param group character(1), name of the column of \code{colData(sce)}
#'   with the replicate of each cell, such as the patient; the replicates
#'   are the folds of the leave-one-out cross-validation
#' @param p_cutoff Cut-off for statistical significance. A horizontal line
#'                 will be drawn at -log10(p_cutoff). Defaults to 0.01
#'                 divided by the number of feature groups. With few
#'                 patients, permutation p-values cannot fall below this
#'                 default, because the smallest attainable p-value is one
#'                 over the number of relabelings.
#' @param fc_cutoff Cut-off for absolute log2 fold-change. A vertical lines
#'                  will be drawn at fc_cutoff.
#' @param assay_type character(1), name of the assay with the features
#' @param stats Output of \code{\link{calculateStats}}; if `NULL`, it is
#'              computed from `sce`
#' @param ... Further arguments to \code{\link{calculateStats}}, such as
#'            `method` or `n_permutations`
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
#' # p-values from the t-distribution keep the example fast; the default
#' # permutation test refits all models for every relabeling of the patients
#' volcanoPlot(sce_single,
#'     target = "Drug", group = "Patient",
#'     p_cutoff = 0.05, fc_cutoff = 0.5, method = "t"
#' )
#'
volcanoPlot <- function(sce,
                        target,
                        group,
                        p_cutoff = NULL,
                        fc_cutoff = 1.0,
                        assay_type = "tfmfeatures",
                        stats = NULL,
                        ...) {
    check_prevalidated(sce, target, group)
    check_number(fc_cutoff, "fc_cutoff")
    if (!is.null(p_cutoff)) check_number(p_cutoff, "p_cutoff", min = 0)
    if (is.null(stats)) {
        stats <- calculateStats(sce, target, group, assay_type, ...)
    }

    if (is.null(p_cutoff)) {
        p_cutoff <- 0.01 / nrow(stats)
    }

    stats |>
        mutate(Feature = ifelse(pvalue < p_cutoff & log2FoldChange > fc_cutoff,
            Feature, ""
        )) |>
        ggplot(aes(log2FoldChange, -log10(pvalue), label = Feature)) +
        geom_vline(
            xintercept = c(0, fc_cutoff), alpha = 0.5,
            linetype = "dashed"
        ) +
        geom_hline(
            yintercept = c(0, -log10(p_cutoff)), alpha = 0.5,
            linetype = "dashed"
        ) +
        xlab("log2 fold change") +
        ylab("-log10 p-value") +
        geom_point() +
        geom_text_repel()
}

#' Plot ROC curves
#'
#' @importFrom ggplot2 ggplot geom_path geom_abline coord_equal
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous xlab ylab ggtitle
#' @importFrom ggplot2 facet_wrap
#' @importFrom dplyr group_by mutate pull
#' @importFrom yardstick roc_auc roc_curve
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param group character(1), name of the column of \code{colData(sce)}
#'   with the replicate of each cell, such as the patient; the replicates
#'   are the folds of the leave-one-out cross-validation
#' @param assay_type character(1), name of the assay with the features
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
#' plotROC(sce_single, target = "Drug", group = "Patient")
#'
plotROC <- function(sce,
                    target,
                    group,
                    assay_type = "tfmfeatures") {
    check_prevalidated(sce, target, group)
    check_assay(sce, assay_type)
    meta_vars <- c(target, group)
    y_hat <- aggregateYhat(sce, target, group, assay_type)
    result <- y_hat |>
        pivot_longer(
            cols = -c(all_of(meta_vars)),
            names_to = "features", values_to = "pred"
        ) |>
        droplevels()
    result[[target]] <- factor(result[[target]], levels = levels(sce[[target]]))

    # title
    reference_level <- levels(result[[target]])[1]
    interest_level <- levels(result[[target]])[2]
    mean_auc <- result |>
        group_by(features) |>
        yardstick::roc_auc(all_of(target), "pred", event_level = "second") |>
        mutate(.estimate = round(.estimate, digits = 3)) |>
        pull(.estimate) |>
        mean() |>
        round(digits = 3)
    title_str <- paste0(
        "Predict ", interest_level, " vs ", reference_level,
        " (mean AUC = ", mean_auc, ")"
    )

    result |>
        group_by(features) |>
        yardstick::roc_curve(all_of(target), "pred", event_level = "second") |>
        ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        geom_path() +
        geom_abline(lty = 3) +
        coord_equal() +
        scale_x_continuous(breaks = c(0, 0.5, 1)) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        xlab("false positive rate (1 - specificity)") +
        ylab("true positive rate (sensitivity)") +
        ggtitle(title_str) +
        facet_wrap(~features)
}

#' Plot AUC comparison
#'
#' @importFrom stats median
#' @importFrom dplyr desc
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot xlab ggtitle geom_jitter aes
#' @importFrom dplyr group_by bind_rows summarize arrange pull
#' @importFrom yardstick roc_auc
#' @importFrom stringr str_remove
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param target character(1), name of the column of \code{colData(sce)}
#'   with the condition of each cell, such as the drug
#' @param group character(1), name of the column of \code{colData(sce)}
#'   with the replicate of each cell, such as the patient; the replicates
#'   are the folds of the leave-one-out cross-validation
#' @param assay_type character(1), name of the assay with the features
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' cell_file <- generateData()
#' sce <- loadData(cell_file)
#' sce <- transformScale(sce)
#'
#' sce$Drug <- as.factor(sce$Drug)
#' sce$Drug <- relevel(sce$Drug, ref = "D1")
#' types <- c("AreaShape", "Intensity", "Texture")
#'
#' sce_single <- predictLOO(
#'     sce,
#'     target = "Drug", group = "Patient",
#'     interest_level = "D7", reference_level = "D1",
#'     types = types,
#'     n_threads = 1
#' )
#'
#' plotAUC(sce_single, target = "Drug", group = "Patient")
#'
plotAUC <- function(sce,
                    target,
                    group,
                    assay_type = "tfmfeatures") {
    check_prevalidated(sce, target, group)
    check_assay(sce, assay_type)
    meta_vars <- c(target, group)
    y_hat <- aggregateYhat(sce, target, group, assay_type)
    result <- y_hat |>
        pivot_longer(
            cols = -c(all_of(meta_vars)),
            names_to = "features", values_to = "pred"
        ) |>
        droplevels()
    result[[target]] <- factor(result[[target]], levels = levels(sce[[target]]))

    reference_level <- levels(result[[target]])[1]
    interest_level <- levels(result[[target]])[2]
    title_str <- paste0("Predict ", interest_level, " vs ", reference_level)

    aucs <- result |>
        group_by(features) |>
        yardstick::roc_auc(all_of(target), "pred", event_level = "second")

    fct_order <- aucs |>
        group_by(features) |>
        summarize(median_auc = median(.estimate)) |>
        arrange(desc(median_auc)) |>
        pull(features)

    aucs |>
        mutate(features = factor(features, levels = fct_order)) |>
        ggplot(aes(.estimate, features)) +
        geom_point() +
        xlab("AUC") +
        ggtitle(title_str)
}
