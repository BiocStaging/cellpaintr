# Provenance of the files in inst/extdata

## header.csv

The header row, without any measurements, of a per-cell table that
CellProfiler's ExportToSpreadsheet module wrote for a live-cell functional
panel of patient-derived cells in our lab (Center of Experimental
Rheumatology, University Hospital Zurich). It lists the metadata columns
(`ImageNumber`, `ObjectNumber`, `Metadata_Patient`, `Metadata_PatientID`,
`Metadata_Well`) and 770 measurement names of the Cells compartment, for four
illumination-corrected channels (prefix `Corr`): `Blue` (nuclear stain),
`CytoID` (CYTO-ID autophagy dye), `Lyso` (LysoTracker), and `TMRM`
(mitochondrial membrane potential).
<!-- TODO(authors): confirm the stains and add the CellProfiler version and pipeline. -->

`generateData()` samples feature names from this header, so that the
simulated table has realistic CellProfiler names and the transformations of
`transformScale()` apply as they would to real data. The file contains no
patient data.

## jump_subset.csv.gz

A small subset of the JUMP Cell Painting dataset (cpg0016) from the Cell
Painting Gallery, released under CC0: 240 U2OS cells from four compound
plates of source 15, treated with DMSO, AMG900, or dexamethasone, with all
CellProfiler measurements of the Cells compartment. The script
`make-jump-subset.R` in this folder created it and documents the selection
and every change to the CellProfiler export.
