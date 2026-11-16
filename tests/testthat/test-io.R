test_that("generateData writes a CellProfiler-style table", {
    set.seed(23)
    cell_file <- generateData()
    expect_true(file.exists(cell_file))

    data <- readr::read_csv(cell_file, show_col_types = FALSE)
    expect_equal(nrow(data), 1000)
    expect_true(all(c(
        "ObjectNumber", "ImageNumber", "Metadata_Well",
        "Metadata_Patient", "Metadata_Drug"
    ) %in% names(data)))
    # five metadata columns and 10 features of each of three types
    expect_equal(ncol(data), 35)
    expect_equal(sum(startsWith(names(data), "AreaShape_")), 10)
    expect_setequal(unique(data$Metadata_Patient), paste0("P", 1:6))
    expect_setequal(unique(data$Metadata_Drug), paste0("D", 1:8))
})

test_that("loadData splits metadata and features", {
    set.seed(23)
    cell_file <- generateData()
    data <- readr::read_csv(cell_file, show_col_types = FALSE)
    sce <- loadData(cell_file)

    expect_s4_class(sce, "SingleCellExperiment")
    expect_equal(dim(sce), c(30L, 1000L))
    expect_setequal(
        names(SummarizedExperiment::colData(sce)),
        c(
            "ObjectNumber", "ImageNumber", "Well", "Patient", "Drug",
            "Location_Center_X", "Location_Center_Y", "Location_Center_Z"
        )
    )
    # features keep their values, one column per cell
    features <- SummarizedExperiment::assay(sce, "features")
    expect_equal(unname(features[, 1]), unname(unlist(data[1, rownames(sce)])))
    expect_equal(sce$Patient, data$Metadata_Patient)
    expect_true(all(is.na(sce$Location_Center_X)))
})

test_that("loadData reads the JUMP subset", {
    cell_file <- system.file("extdata", "jump_subset.csv.gz",
        package = "cellpaintr"
    )
    skip_if(cell_file == "")
    sce <- loadData(cell_file)
    expect_equal(ncol(sce), 240L)
    expect_setequal(
        unique(sce$Compound), c("DMSO", "AMG900", "dexamethasone")
    )
    expect_equal(length(unique(sce$Plate)), 4L)
    expect_true(all(c("Well", "Site") %in% names(colData(sce))))
    expect_true(nrow(sce) > 1000)
})

test_that("loadData rejects tables that are not CellProfiler exports", {
    expect_error(loadData("no_such_file.csv"), "does not exist")
    expect_error(loadData(c("a.csv", "b.csv")), "single non-empty string")

    write_table <- function(data) {
        file <- tempfile(fileext = ".csv")
        readr::write_csv(data, file)
        file
    }
    expect_error(
        loadData(write_table(data.frame(
            Metadata_Well = "A01", AreaShape_Area = 1
        ))),
        "ImageNumber"
    )
    expect_error(
        loadData(write_table(data.frame(
            ImageNumber = 1, AreaShape_Area = 1
        ))),
        "Metadata_"
    )
    expect_error(
        loadData(write_table(data.frame(
            ImageNumber = 1, Metadata_Well = "A01", Drug = "D1",
            AreaShape_Area = 1
        ))),
        "numeric"
    )
})
