test_that("remove missing values", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)

    assay(sce)[1, 1] <- NA
    sce <- removeNAs(sce)
    n_nas <- sum(is.na(assay(sce)))

    expect_equal(n_nas, 0)
})

test_that("remove outliers", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)

    sce <- removeOutliers(sce, min = 0, max = 17)

    expect_equal(sum(table(sce$ImageNumber) > 17), 0)
})

test_that("remove low variance features (robust = FALSE)", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)

    assay(sce)["AreaShape_Zernike_2_2", ] <- 0
    sce <- removeLowVariance(sce, robust = FALSE)

    expect_disjoint("AreaShape_Zernike_2_2", rownames(sce))
})

test_that("remove low variance features (robust = TRUE)", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)

    assay(sce)["AreaShape_Zernike_2_2", ] <- 0
    sce <- removeLowVariance(sce, robust = TRUE)

    expect_disjoint("AreaShape_Zernike_2_2", rownames(sce))
})

test_that("remove zero inflation features", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)

    assay(sce)["AreaShape_Zernike_2_2", ] <- 0
    sce <- removeZeroInflation(sce)

    expect_disjoint("AreaShape_Zernike_2_2", rownames(sce))
})

test_that("remove row and col effects", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)

    # check input
    expect_error(
        removeRowColEffect(sce),
        "'Col' in `colData\\(sce\\)\\$Col` must be a factor"
    )
    col <- substr(sce$Well, 2, 3)
    col <- factor(as.numeric(col))
    sce$Col <- col
    expect_error(
        removeRowColEffect(sce),
        "'Row' in `colData\\(sce\\)\\$Row` must be a factor"
    )

    # check run
    row <- substr(sce$Well, 1, 1)
    row <- factor(row, levels = sort(unique(row)))
    sce$Row <- row
    sce_corrected <- removeRowColEffect(sce)
    expect_s4_class(sce_corrected, "SingleCellExperiment")
    expect_true("corrected" %in% assayNames(sce_corrected))
    expect_equal(dim(sce_corrected), dim(sce))
})

test_that("filters use the assay given by assay_type", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    assay(sce, "copy") <- assay(sce, "features")
    assay(sce, "copy")["AreaShape_Zernike_2_2", ] <- 0
    assay(sce, "copy")[1, 1] <- NA

    expect_true("AreaShape_Zernike_2_2" %in%
        rownames(removeLowVariance(sce)))
    expect_disjoint(
        "AreaShape_Zernike_2_2",
        rownames(removeZeroInflation(sce, assay_type = "copy"))
    )
    expect_equal(ncol(removeNAs(sce)), ncol(sce))
    expect_equal(ncol(removeNAs(sce, assay_type = "copy")), ncol(sce) - 1)
    expect_error(removeNAs(sce, assay_type = "nope"), "not found")
})

test_that("transformScale standardizes and caps infinite values", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)

    standard <- assay(transformScale(sce), "tfmfeatures")
    expect_equal(unname(rowMeans(standard)), rep(0, nrow(sce)),
        tolerance = 1e-10
    )
    expect_equal(unname(apply(standard, 1, sd)), rep(1, nrow(sce)),
        tolerance = 1e-10
    )

    robust <- assay(transformScale(sce, robust = TRUE), "tfmfeatures")
    expect_equal(unname(apply(robust, 1, median)), rep(0, nrow(sce)),
        tolerance = 1e-10
    )
    expect_equal(unname(apply(robust, 1, mad)), rep(1, nrow(sce)),
        tolerance = 1e-10
    )

    # intensities are log-transformed, so a zero becomes -Inf and is capped
    # at the smallest finite value of its feature
    intensity <- grep("^Intensity_MeanIntensity_", rownames(sce),
        value = TRUE
    )[1]
    assay(sce)[intensity, 1] <- 0
    capped <- assay(transformScale(sce), "tfmfeatures")
    expect_true(all(is.finite(capped)))
    expect_equal(unname(capped[intensity, 1]), min(capped[intensity, -1]))

    # the transformation is monotone within a feature
    original <- assay(sce)[intensity, -1]
    expect_equal(
        order(capped[intensity, -1]), order(original),
        ignore_attr = TRUE
    )
})

test_that("preprocessing functions validate their arguments", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)

    expect_error(removeNAs(assay(sce)), "SingleCellExperiment")
    expect_error(removeLowVariance(sce, robust = "yes"), "TRUE")
    expect_error(removeZeroInflation(sce, proportion = -1), "at least 0")
    expect_error(removeOutliers(sce, min = 10, max = 5), "at least 10")
    expect_error(transformScale(sce, assay_type = "tfmfeatures"), "not found")
    expect_error(plotCellsPerImage(sce, bins = 0), "at least 1")
    expect_error(plotPCACor(transformScale(sce)), "runPCA")
})

test_that("removeRowColEffect leaves out empty wells", {
    # an additive plate: every cell is the sum of a row and a column effect,
    # with well A1 empty
    grid <- expand.grid(row = LETTERS[1:4], col = 1:4, cell = 1:3)
    grid <- grid[!(grid$row == "A" & grid$col == 1), ]
    row_effect <- c(A = 1, B = 2, C = 4, D = 8)
    col_effect <- c(0, 10, 20, 30)
    value <- row_effect[as.character(grid$row)] + col_effect[grid$col]
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(tfmfeatures = rbind(feature = value))
    )
    sce$Row <- factor(grid$row)
    sce$Col <- factor(grid$col)

    # the median polish recovers the effects exactly, so nothing is left;
    # setting the empty well to zero would leave residuals
    corrected <- assay(removeRowColEffect(sce), "corrected")
    expect_equal(unname(corrected[1, ]), rep(0, ncol(sce)), tolerance = 1e-8)
})
