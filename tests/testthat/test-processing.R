test_that("train prediction model on types", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)

    sce$Drug <- as.factor(sce$Drug)
    sce$Drug <- relevel(sce$Drug, ref = "D1")
    types <- c("AreaShape", "Intensity", "Texture")

    # leave-one-out score
    sce_single <- predictLOO(
        sce,
        target = "Drug", group = "Patient",
        interest_level = "D7", reference_level = "D1",
        types = types,
        n_threads = 1
    )
    scores <- reducedDim(sce_single, "prevalidated")

    expect_setequal(colnames(scores), c("all", types))
    expect_gte(min(scores), 0)
    expect_lte(max(scores), 1)
})

test_that("train prediction model on types with weights", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)

    # aggregate
    aggr <- scrapper::aggregateAcrossCells(
        assay(sce, "tfmfeatures"),
        factors = colData(sce)[, c("ImageNumber", "Drug", "Patient")]
    )
    sce_aggr <- SingleCellExperiment(
        assays = list(sums = aggr$sums),
        colData = DataFrame(aggr$combinations, ncells = aggr$counts)
    )

    sce_aggr$Drug <- as.factor(sce_aggr$Drug)
    sce_aggr$Drug <- relevel(sce_aggr$Drug, ref = "D1")
    types <- c("AreaShape", "Intensity", "Texture")

    # leave-one-out score
    sce_single <- predictLOO(
        sce_aggr,
        weights = "ncells",
        target = "Drug", group = "Patient",
        interest_level = "D7", reference_level = "D1",
        types = types,
        n_threads = 1,
        assay_type = "sums"
    )
    scores <- reducedDim(sce_single, "prevalidated")

    expect_setequal(colnames(scores), c("all", types))
    expect_gte(min(scores), 0)
    expect_lte(max(scores), 1)
})

test_that("train prediction model on channels", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)

    sce$Drug <- as.factor(sce$Drug)
    sce$Drug <- relevel(sce$Drug, ref = "D1")
    channels <- c("CorrBlue", "CorrTMRM", "CorrLyso", "CorrCytoID")

    # leave-one-out score
    sce_single <- predictLOO(
        sce,
        target = "Drug", group = "Patient",
        interest_level = "D7", reference_level = "D1",
        channels = channels,
        n_threads = 1
    )
    scores <- reducedDim(sce_single, "prevalidated")

    expect_setequal(colnames(scores), c("all", channels))
    expect_gte(min(scores), 0)
    expect_lte(max(scores), 1)
})

test_that("extract prediction scores", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)

    sce$Drug <- as.factor(sce$Drug)
    sce$Drug <- relevel(sce$Drug, ref = "D1")
    types <- c("AreaShape", "Intensity", "Texture")

    # leave-one-out score
    sce_single <- predictLOO(
        sce,
        target = "Drug", group = "Patient",
        interest_level = "D7", reference_level = "D1",
        types = types,
        n_threads = 1
    )
    stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    )

    expect_setequal(
        colnames(stats),
        c("Target", "Feature", "pvalue", "log2FoldChange")
    )
    expect_equal(nrow(stats), 4)
    expect_equal(ncol(stats), 4)
    expect_gte(min(stats$pvalue), 0)
    expect_lte(max(stats$pvalue), 1)
})

test_that("predictLOO validates its arguments", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)
    fit <- function(...) {
        args <- list(
            sce = sce, target = "Drug", group = "Patient",
            interest_level = "D7", reference_level = "D1",
            progress = FALSE
        )
        args[names(list(...))] <- list(...)
        do.call(predictLOO, args)
    }

    expect_error(fit(target = "Dose"), "not a column")
    expect_error(fit(group = "Drug"), "different columns")
    expect_error(fit(interest_level = "D9"), "Not found")
    expect_error(fit(reference_level = "D7"), "must differ")
    expect_error(fit(types = 1), "character vector")
    expect_error(fit(weights = "w"), "not a column")
    expect_error(fit(n_threads = 0), "at least 1")
    expect_error(fit(assay_type = "corrected"), "not found")
    single <- sce[, sce$Patient == "P1"]
    expect_error(fit(sce = single), "at least two levels")
})

test_that("predictLOO keeps cell names and accepts a character target", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)
    colnames(sce) <- paste0("cell", seq_len(ncol(sce)))
    sce <- sce[, sce$Patient %in% c("P1", "P2", "P3")]
    expect_type(sce$Drug, "character")

    sce_single <- predictLOO(sce,
        target = "Drug", group = "Patient",
        interest_level = "D7", reference_level = "D1",
        progress = FALSE
    )
    expect_s3_class(sce_single$Drug, "factor")
    expect_equal(
        colnames(sce_single),
        colnames(sce)[sce$Drug %in% c("D1", "D7")]
    )
    expect_equal(
        rownames(reducedDim(sce_single, "prevalidated")),
        colnames(sce_single)
    )
})

test_that("functions on predictLOO output check for its scores", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    expect_error(
        calculateStats(sce, target = "Drug", group = "Patient"),
        "prevalidated"
    )
    expect_error(plotLOO(sce, target = "Drug", group = "Patient"),
        "prevalidated"
    )
})
