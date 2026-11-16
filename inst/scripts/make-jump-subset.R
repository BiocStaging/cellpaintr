# Provenance of inst/extdata/jump_subset.csv.gz
#
# A small subset of the JUMP Cell Painting dataset (cpg0016; Chandrasekaran et
# al. 2023, doi:10.1101/2023.03.23.534023) from the Cell Painting Gallery
# (Weisbart et al. 2024, doi:10.1038/s41592-024-02399-z), released under CC0.
# U2OS cells were stained with the Cell Painting protocol (Hoechst for DNA,
# concanavalin A for the ER, SYTO 14 for RNA, phalloidin and wheat germ
# agglutinin for actin, Golgi, and plasma membrane (AGP), and MitoTracker for
# mitochondria), imaged, and measured with CellProfiler; this subset keeps the
# per-cell measurements of the Cells compartment as CellProfiler exported
# them (Cells.csv), joined with image metadata from Image.csv and compound
# names from the JUMP metadata.
#
# Selection: the first four compound plates of batch 2021_12_03_Batch1 of
# source 15, the DMSO negative control and two of the eight JUMP positive
# controls, AMG900 (a strong effect) and dexamethasone (a weak effect). For
# each plate and compound, two random wells; for each well, one random site
# and 10 random cells of that site. Plates take the place of patients: every
# plate carries all three compounds, so the design is paired.
#
# Each Cells.csv holds all cells of a plate (about 27 GB). Its rows are sorted
# by ImageNumber, so the script finds the rows of each selected site by binary
# search with HTTP range requests and downloads only those, a few megabytes
# per site. Image.csv (about 115 MB per plate) is downloaded in full.
#
# Changes to the CellProfiler export:
#   - Metadata_Plate, Metadata_Well, Metadata_Site (from Image.csv) and
#     Metadata_Compound (from the JUMP metadata) are added, as CellProfiler's
#     ExportToSpreadsheet module does with "Add image metadata columns to your
#     object data file"
#   - ImageNumber is renumbered so that it is unique across plates
#   - numeric values are rounded to four significant digits to keep the file
#     small
#
# Usage: Rscript make-jump-subset.R [output file]
# Needs readr, dplyr, and curl.

suppressPackageStartupMessages({
    library(dplyr)
})

out_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_file)) out_file <- "jump_subset.csv.gz"
cache_dir <- file.path(tempdir(), "jump_subset")
dir.create(cache_dir, showWarnings = FALSE)
set.seed(2026)

source_id <- "source_15"
batch_id <- "2021_12_03_Batch1"
n_plates <- 4
compounds <- c("DMSO", "AMG900", "dexamethasone")
wells_per_compound <- 2
cells_per_well <- 10

metadata_url <- paste0(
    "https://raw.githubusercontent.com/jump-cellpainting/datasets/",
    "016e865fa0691244e0860943e41c7d6a88ed2580/metadata/"
)
analysis_url <- paste0(
    "https://cellpainting-gallery.s3.amazonaws.com/cpg0016-jump/",
    source_id, "/workspace/analysis/", batch_id, "/"
)

# 1. Plates and wells ---------------------------------------------------------

read_metadata <- function(file) {
    readr::read_csv(paste0(metadata_url, file), show_col_types = FALSE)
}
plates <- read_metadata("plate.csv.gz")
wells <- read_metadata("well.csv.gz")
controls <- read_metadata("perturbation_control.csv")

layout <- wells |>
    inner_join(plates, by = c("Metadata_Source", "Metadata_Plate")) |>
    filter(
        Metadata_Source == source_id, Metadata_Batch == batch_id,
        Metadata_PlateType == "COMPOUND"
    ) |>
    inner_join(select(controls, Metadata_JCP2022, Metadata_Name),
        by = "Metadata_JCP2022"
    ) |>
    filter(Metadata_Name %in% compounds)

selected_plates <- layout |>
    count(Metadata_Plate, Metadata_Name) |>
    tidyr::pivot_wider(
        names_from = Metadata_Name, values_from = n, values_fill = 0
    ) |>
    filter(if_all(all_of(compounds), ~ .x >= wells_per_compound)) |>
    arrange(Metadata_Plate) |>
    slice_head(n = n_plates) |>
    pull(Metadata_Plate)
stopifnot(length(selected_plates) == n_plates)

selected_wells <- layout |>
    filter(Metadata_Plate %in% selected_plates) |>
    group_by(Metadata_Plate, Metadata_Name) |>
    slice_sample(n = wells_per_compound) |>
    ungroup() |>
    arrange(Metadata_Plate, Metadata_Well)

# 2. Byte ranges of a sorted CSV ----------------------------------------------

fetch <- function(url, from, to) {
    handle <- curl::new_handle()
    range <- sprintf("bytes=%.0f-%.0f", from, to)
    curl::handle_setheaders(handle, Range = range)
    response <- curl::curl_fetch_memory(url, handle = handle)
    if (!response$status_code %in% c(200, 206)) {
        stop("HTTP ", response$status_code, " for ", url)
    }
    rawToChar(response$content)
}

remote_size <- function(url) {
    handle <- curl::new_handle(nobody = TRUE)
    response <- curl::curl_fetch_memory(url, handle = handle)
    headers <- curl::parse_headers_list(response$headers)
    as.numeric(headers[["content-length"]])
}

chunk <- 60000 # bytes; a row of Cells.csv has about 37,000

# ImageNumber of the first complete row at or after byte `offset`
image_at <- function(url, offset) {
    text <- fetch(url, offset, offset + chunk)
    start <- regexpr("\n", text)
    rest <- substring(text, start + 1)
    as.integer(sub(",.*", "", rest))
}

# first byte offset whose next complete row has ImageNumber >= `image`
find_offset <- function(url, size, image) {
    lo <- 0
    hi <- size
    while (hi - lo > chunk) {
        mid <- floor((lo + hi) / 2)
        if (image_at(url, mid) < image) lo <- mid else hi <- mid
    }
    lo
}

# all rows of one image, read forward from `offset` in blocks
read_image_rows <- function(url, size, offset, image, header) {
    block <- 4e6
    text <- ""
    position <- offset
    repeat {
        last <- min(position + block, size) - 1
        text <- paste0(text, fetch(url, position, last))
        position <- position + block
        lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
        # drop the partial first line and the possibly partial last line
        lines <- lines[-1]
        complete <- if (position >= size) lines else head(lines, -1)
        images <- as.integer(sub(",.*", "", complete))
        if (any(images > image, na.rm = TRUE) || position >= size) break
    }
    rows <- complete[!is.na(images) & images == image]
    cells <- readr::read_csv(I(c(header, rows)), show_col_types = FALSE,
        col_types = readr::cols(.default = readr::col_double())
    )
    # every row must have one value per column of the header
    stopifnot(nrow(readr::problems(cells)) == 0)
    cells
}

# 3. Cells of the selected sites ----------------------------------------------

cells <- lapply(selected_plates, function(plate) {
    base <- paste0(analysis_url, plate, "/analysis/", plate, "/")
    image_file <- file.path(cache_dir, paste0(plate, "_Image.csv"))
    if (!file.exists(image_file)) {
        curl::curl_download(paste0(base, "Image.csv"), image_file)
    }
    images <- readr::read_csv(image_file,
        col_select = c("ImageNumber", "Metadata_Well", "Metadata_Site"),
        show_col_types = FALSE
    )
    cells_url <- paste0(base, "Cells.csv")
    size <- remote_size(cells_url)
    # the header line alone is longer than one chunk
    header_text <- fetch(cells_url, 0, 10 * chunk)
    stopifnot(grepl("\n", header_text))
    header <- sub("\n.*", "", header_text)

    plate_wells <- filter(selected_wells, Metadata_Plate == plate)
    lapply(seq_len(nrow(plate_wells)), function(i) {
        well <- plate_wells$Metadata_Well[i]
        candidates <- sample(images$ImageNumber[images$Metadata_Well == well])
        for (image in candidates) {
            offset <- find_offset(cells_url, size, image)
            rows <- read_image_rows(cells_url, size, offset, image, header)
            if (nrow(rows) >= cells_per_well) break
        }
        stopifnot(nrow(rows) >= cells_per_well)
        message(plate, " ", well, ": image ", image, ", ", nrow(rows), " cells")
        rows |>
            slice_sample(n = cells_per_well) |>
            mutate(
                Metadata_Plate = plate,
                Metadata_Well = well,
                Metadata_Site =
                    images$Metadata_Site[images$ImageNumber == image],
                Metadata_Compound = plate_wells$Metadata_Name[i],
                .after = ObjectNumber
            )
    }) |> bind_rows()
}) |> bind_rows()

# 4. Output -------------------------------------------------------------------

cells <- cells |>
    mutate(
        ImageNumber = match(
            paste(Metadata_Plate, ImageNumber),
            unique(paste(Metadata_Plate, ImageNumber))
        ),
        across(
            where(is.double) & !c(ImageNumber, ObjectNumber),
            ~ signif(.x, 4)
        )
    )
readr::write_csv(cells, out_file)
message(
    "Wrote ", nrow(cells), " cells and ", ncol(cells), " columns to ", out_file,
    " (", round(file.size(out_file) / 1e6, 2), " MB)"
)
