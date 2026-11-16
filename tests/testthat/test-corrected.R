simulate_paired <- function() {
    set.seed(23)
    cell_file <- generate_data()
    sce <- loadData(cell_file)
    sce <- transformScale(sce)
    sce$Drug <- relevel(as.factor(sce$Drug), ref = "D1")
    sce
}

fit_paired <- function(sce, patients) {
    set.seed(1)
    predictLOO(sce[, sce$Patient %in% patients],
        target = "Drug", group = "Patient",
        interest_level = "D8", reference_level = "D1",
        types = "Texture", n_threads = 1, progress = FALSE
    )
}

test_that("correction factor follows the pairwise-interaction model", {
    expect_equal(correction_factor(4), sqrt(3))
    expect_equal(correction_factor(6), sqrt(10 / 4))
    # equicorrelated score differences with correlation 1 / (n - 1)
    n <- 7
    rho <- 1 / (n - 1)
    expect_equal(correction_factor(n)^2, (1 + (n - 1) * rho) / (1 - rho))
})

test_that("corrected t-test divides the paired t-statistic by the factor", {
    sce <- simulate_paired()
    sce_single <- fit_paired(sce, c("P1", "P2", "P3", "P4"))
    stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "corrected"
    )
    t_stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    )
    statistic <- qt(t_stats$pvalue, 3, lower.tail = FALSE)

    expect_equal(stats$Feature, c("all", "Texture"))
    expect_equal(stats$log2FoldChange, t_stats$log2FoldChange)
    expect_equal(
        stats$pvalue,
        pt(statistic / sqrt(3), 3, lower.tail = FALSE)
    )
    # the correction makes the test more conservative in the upper tail
    expect_true(all(stats$pvalue[statistic > 0] > t_stats$pvalue[statistic > 0]))
})

test_that("corrected t-test needs no settings of predictLOO()", {
    sce <- simulate_paired()
    sce_single <- fit_paired(sce, c("P1", "P2", "P3"))
    S4Vectors::metadata(sce_single)$predictLOO <- NULL
    expect_no_error(calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "corrected"
    ))
})

test_that("corrected t-test needs a paired design", {
    sce <- simulate_paired()
    keep <- (sce$Patient %in% c("P1", "P2", "P3") & sce$Drug == "D1") |
        (sce$Patient %in% c("P4", "P5", "P6") & sce$Drug == "D8")
    set.seed(1)
    sce_single <- predictLOO(sce[, keep],
        target = "Drug", group = "Patient",
        interest_level = "D8", reference_level = "D1",
        n_threads = 1, progress = FALSE
    )
    expect_error(
        calculateStats(sce_single,
            target = "Drug", group = "Patient", method = "corrected"
        ),
        "paired design"
    )
})

test_that("Nadeau-Bengio divides the paired t-statistic by its factor", {
    expect_equal(correction_factor(4, "nadeau_bengio"), sqrt(7 / 3))
    # equicorrelated score differences with correlation 1 / n
    n <- 7
    rho <- 1 / n
    expect_equal(
        correction_factor(n, "nadeau_bengio")^2,
        (1 + (n - 1) * rho) / (1 - rho)
    )
    # it corrects less than the pairwise-interaction factor
    expect_true(all(
        correction_factor(3:20, "nadeau_bengio") < correction_factor(3:20)
    ))

    sce <- simulate_paired()
    sce_single <- fit_paired(sce, c("P1", "P2", "P3", "P4"))
    stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "nadeau_bengio"
    )
    t_stats <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "t"
    )
    corrected <- calculateStats(sce_single,
        target = "Drug", group = "Patient", method = "corrected"
    )
    statistic <- qt(t_stats$pvalue, 3, lower.tail = FALSE)
    expect_equal(
        stats$pvalue,
        pt(statistic / sqrt(7 / 3), 3, lower.tail = FALSE)
    )
    positive <- statistic > 0
    expect_true(all(stats$pvalue[positive] > t_stats$pvalue[positive]))
    expect_true(all(stats$pvalue[positive] < corrected$pvalue[positive]))
})

test_that("Nadeau-Bengio needs a paired design and works with two patients", {
    sce <- simulate_paired()
    expect_no_error(calculateStats(fit_paired(sce, c("P1", "P2")),
        target = "Drug", group = "Patient", method = "nadeau_bengio"
    ))
    keep <- (sce$Patient %in% c("P1", "P2", "P3") & sce$Drug == "D1") |
        (sce$Patient %in% c("P4", "P5", "P6") & sce$Drug == "D8")
    set.seed(1)
    sce_single <- predictLOO(sce[, keep],
        target = "Drug", group = "Patient",
        interest_level = "D8", reference_level = "D1",
        n_threads = 1, progress = FALSE
    )
    expect_error(
        calculateStats(sce_single,
            target = "Drug", group = "Patient", method = "nadeau_bengio"
        ),
        "Nadeau-Bengio t-test needs a paired design"
    )
})

test_that("corrected t-test needs at least three patients", {
    sce <- simulate_paired()
    sce_single <- fit_paired(sce, c("P1", "P2"))
    expect_error(
        calculateStats(sce_single,
            target = "Drug", group = "Patient", method = "corrected"
        ),
        "at least 3"
    )
})
