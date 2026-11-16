# small designs keep the number of refits low
simulate_sce <- function() {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)
    sce$Drug <- relevel(as.factor(sce$Drug), ref = "D1")
    sce
}

fit_loo <- function(sce) {
    predictLOO(sce,
        target = "Drug", group = "Patient",
        interest_level = "D8", reference_level = "D1",
        n_threads = 1, progress = FALSE
    )
}

on_grid <- function(pvalue, n) all(abs(pvalue * n - round(pvalue * n)) < 1e-8)

test_that("predictLOO records its settings", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3")])
    settings <- S4Vectors::metadata(sce_single)$predictLOO

    expect_equal(settings$target, "Drug")
    expect_equal(settings$group, "Patient")
    expect_equal(settings$interest_level, "D8")
    expect_equal(settings$reference_level, "D1")
    expect_equal(settings$assay_type, "tfmfeatures")
})

test_that("permutation test enumerates all swaps in a paired design", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3")])
    expect_message(
        stats <- calculateStats(sce_single, target = "Drug", group = "Patient"),
        "all 4 relabelings"
    )

    expect_setequal(
        colnames(stats),
        c("Target", "Feature", "pvalue", "log2FoldChange")
    )
    # P1 keeps its labels, so 2^(3 - 1) relabelings
    expect_equal(attr(stats, "relabelings"), 4)
    expect_true(on_grid(stats$pvalue, 4))
    expect_gte(min(stats$pvalue), 1 / 4)
    expect_lte(max(stats$pvalue), 1)
})

test_that("permutation test samples swaps when there are too many", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3")])
    stats <- suppressMessages(calculateStats(sce_single,
        target = "Drug", group = "Patient", n_permutations = 2
    ))

    expect_equal(attr(stats, "relabelings"), 3)
    expect_true(on_grid(stats$pvalue, 3))
})

test_that("permutation test permutes patients in an unpaired design", {
    sce <- simulate_sce()
    keep <- (sce$Patient %in% c("P1", "P2") & sce$Drug == "D1") |
        (sce$Patient %in% c("P3", "P4") & sce$Drug == "D8")
    sce_single <- fit_loo(sce[, keep])
    expect_message(
        stats <- calculateStats(sce_single, target = "Drug", group = "Patient"),
        "unpaired design"
    )

    # two levels per side, so P1 keeps its level: choose(3, 1) relabelings
    expect_equal(attr(stats, "relabelings"), choose(3, 1))
    expect_true(on_grid(stats$pvalue, 3))
})

test_that("permutation test handles mixed designs", {
    sce <- simulate_sce()
    # P1 has only D1 cells, P4 only D8 cells, P2 and P3 have both
    keep <- sce$Patient %in% c("P1", "P2", "P3", "P4") &
        !(sce$Patient == "P1" & sce$Drug == "D8") &
        !(sce$Patient == "P4" & sce$Drug == "D1")
    sce_single <- fit_loo(sce[, keep])
    expect_message(
        stats <- calculateStats(sce_single, target = "Drug", group = "Patient"),
        "mixed design"
    )

    # P3 keeps or swaps its levels, P2 keeps them, and either P1 or P4 gets D8
    expect_equal(attr(stats, "relabelings"), 2 * choose(2, 1))
    expect_true(on_grid(stats$pvalue, 4))
})

test_that("permutation test samples relabelings in mixed designs", {
    sce <- simulate_sce()
    keep <- sce$Patient %in% c("P1", "P2", "P3") &
        !(sce$Patient == "P1" & sce$Drug == "D8")
    sce_single <- fit_loo(sce[, keep])
    stats <- suppressMessages(calculateStats(sce_single,
        target = "Drug", group = "Patient", n_permutations = 2
    ))

    expect_equal(attr(stats, "relabelings"), 3)
    expect_true(on_grid(stats$pvalue, 3))
})

test_that("relabelings keep one patient fixed when exchanging all is one", {
    # paired: the first patient never swaps, and the observed labels are left out
    grid <- patient_relabelings("paired", 4, rep(TRUE, 4),
        exact = TRUE, n_permutations = 100, symmetric = TRUE
    )
    expect_equal(nrow(grid), 2^3 - 1)
    expect_false(any(grid[, 1]))
    expect_false(any(duplicated(grid)))
    expect_true(all(rowSums(grid) > 0))
    sampled <- patient_relabelings("paired", 4, rep(TRUE, 4),
        exact = FALSE, n_permutations = 50, symmetric = TRUE
    )
    expect_false(any(sampled[, 1]))

    # without symmetry all 2^n swaps are listed
    grid <- patient_relabelings("paired", 4, rep(TRUE, 4),
        exact = TRUE, n_permutations = 100, symmetric = FALSE
    )
    expect_equal(nrow(grid), 2^4 - 1)

    # unpaired with two patients per level: the first patient keeps its level
    observed <- c(FALSE, TRUE, FALSE, TRUE)
    grid <- patient_relabelings("unpaired", 4, observed,
        exact = TRUE, n_permutations = 100, symmetric = TRUE
    )
    expect_equal(nrow(grid), choose(3, 2) - 1)
    expect_true(all(grid[, 1] == observed[1]))
    expect_true(all(rowSums(grid) == 2))
    sampled <- patient_relabelings("unpaired", 4, observed,
        exact = FALSE, n_permutations = 50, symmetric = TRUE
    )
    expect_true(all(sampled[, 1] == observed[1]))
    expect_true(all(rowSums(sampled) == 2))

    # mixed: the first patient with both levels never swaps
    both <- c(FALSE, TRUE, TRUE, FALSE)
    observed <- c(FALSE, TRUE, TRUE, TRUE)
    grid <- patient_relabelings("mixed", 4, observed,
        exact = TRUE, n_permutations = 100, both = both, symmetric = TRUE
    )
    expect_equal(nrow(grid), 2 * choose(2, 1) - 1)
    expect_false(any(grid[, 2]))
})

test_that("exchanging the levels in all patients keeps the statistic", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3", "P4")])
    settings <- S4Vectors::metadata(sce_single)$predictLOO
    swapped <- sce_single
    swapped$Drug <- factor(
        ifelse(swapped$Drug == "D8", "D1", "D8"),
        levels = levels(swapped$Drug)
    )
    refit <- do.call(
        predictLOO,
        c(list(sce = swapped), settings, list(progress = FALSE))
    )
    statistic <- function(x) {
        score_statistics(x, "Drug", "Patient", settings$assay_type,
            "D8", "D1"
        )$statistic
    }
    refit_same <- do.call(
        predictLOO,
        c(list(sce = sce_single), settings, list(progress = FALSE))
    )

    # the exchange differs from the observed labels by no more than a refit
    # with a new forest seed does
    noise <- abs(statistic(refit_same) - statistic(sce_single))
    expect_lte(
        abs(statistic(refit) - statistic(sce_single)),
        max(3 * noise, 0.1 * abs(statistic(sce_single)))
    )
})

test_that("permutation test needs the settings of predictLOO", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3")])
    S4Vectors::metadata(sce_single)$predictLOO <- NULL

    expect_error(
        calculateStats(sce_single, target = "Drug", group = "Patient"),
        "predictLOO"
    )
    expect_no_error(calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    ))
})

test_that("levels come from predictLOO, not from the factor order", {
    sce <- simulate_sce()
    sce <- sce[, sce$Patient %in% c("P1", "P2", "P3")]
    sce$Drug <- relevel(sce$Drug, ref = "D8")
    sce_single <- fit_loo(sce)
    stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    )

    expect_equal(stats$Target, "D8")
    expect_gt(stats$log2FoldChange, 0)
})

test_that("volcano plot accepts precomputed statistics", {
    sce <- simulate_sce()
    sce_single <- fit_loo(sce[, sce$Patient %in% c("P1", "P2", "P3")])
    stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    )
    p <- volcanoPlot(sce_single,
        target = "Drug", group = "Patient",
        stats = stats, p_cutoff = 0.05
    )

    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_name(p$mapping$x), "log2FoldChange")
})
