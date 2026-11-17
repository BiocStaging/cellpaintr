test_that("plot cells per image", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    p <- plotCellsPerImage(sce)

    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$x), "n")
})

test_that("plot pca cor", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    sce <- transformScale(sce, robust = TRUE)
    sce <- scater::runPCA(sce, exprs_values = "tfmfeatures", ncomponents = 10)
    p <- plotPCACor(sce, filter_by = 1)

    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$x), "PC")
})

test_that("plot prediction scores", {
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
    p <- plotLOO(sce_single, target = "Drug", group = "Patient")
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$y), "value")

    # volcano plot
    p <- volcanoPlot(sce_single,
        target = "Drug", group = "Patient",
        p_cutoff = 0.05, fc_cutoff = 0.5, method = "t"
    )
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$x), "log2FoldChange")

    # roc curve
    p <- plotROC(sce_single, target = "Drug", group = "Patient")
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$y), "sensitivity")

    # auc summary
    p <- plotAUC(sce_single, target = "Drug", group = "Patient")
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$x), ".estimate")

    # no p_cutoff defined
    p <- volcanoPlot(sce_single,
        target = "Drug", group = "Patient",
        fc_cutoff = 0.5, method = "t"
    )
    expect_s3_class(p, "ggplot")
})

test_that("ROC and AUC use the levels of predictLOO, not the factor order", {
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)
    sce <- sce[, sce$Drug %in% c("D1", "D8")]
    # the interest level sorts before the reference level
    sce$Condition <- ifelse(sce$Drug == "D8", "A_drug", "Z_control")

    set.seed(1)
    fit_character <- predictLOO(sce,
        target = "Condition", group = "Patient",
        interest_level = "A_drug", reference_level = "Z_control",
        types = "Texture", progress = FALSE
    )
    expect_equal(levels(fit_character$Condition), c("Z_control", "A_drug"))

    sce$Condition <- factor(sce$Condition, levels = c("A_drug", "Z_control"))
    set.seed(1)
    fit_factor <- predictLOO(sce,
        target = "Condition", group = "Patient",
        interest_level = "A_drug", reference_level = "Z_control",
        types = "Texture", progress = FALSE
    )
    expect_equal(levels(fit_factor$Condition), c("A_drug", "Z_control"))

    auc <- function(fit) {
        plotAUC(fit, target = "Condition", group = "Patient")$data$.estimate
    }
    # D8 shifts the texture features, so its scores separate from the
    # control in the direction of the interest level
    expect_equal(auc(fit_character), auc(fit_factor))
    expect_true(all(auc(fit_factor) > 0.5))
    roc <- plotROC(fit_factor, target = "Condition", group = "Patient")
    expect_match(roc$labels$title, "Predict A_drug vs Z_control")
})
