test_that("permutation test gives the same p-values with a future plan", {
    skip_if_not_installed("future")
    # multisession workers load the installed package, not a dev version
    skip_if_not("cellpaintr" %in% rownames(utils::installed.packages()))
    skip_on_cran()
    set.seed(23)
    cell_file <- generateData()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)
    sce$Drug <- relevel(as.factor(sce$Drug), ref = "D1")
    sce_single <- predictLOO(sce[, sce$Patient %in% c("P1", "P2", "P3")],
        target = "Drug", group = "Patient",
        interest_level = "D8", reference_level = "D1",
        types = "Texture", progress = FALSE
    )

    old_plan <- future::plan(future::sequential)
    on.exit(future::plan(old_plan), add = TRUE)
    set.seed(1)
    sequential <- suppressMessages(calculateStats(sce_single,
        target = "Drug", group = "Patient"
    ))
    future::plan(future::multisession, workers = 2)
    set.seed(1)
    parallel <- suppressMessages(calculateStats(sce_single,
        target = "Drug", group = "Patient"
    ))
    # future.seed gives each relabeling its own reproducible stream
    expect_equal(parallel$pvalue, sequential$pvalue)
})
