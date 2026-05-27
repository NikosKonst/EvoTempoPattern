# EvoTempoPattern

This repository contains all scripts and processed data required to reproduce the single-cell analysis, trajectory inference, temporal transcription factor identification, and figure generation performed in this study, “Evolutionary dynamics of temporal transcription factor series in the insect visual brain”, Filippopoulou et al. 2026.

---

# Repository Structure

## `single_cell_analysis/`

This folder contains one R script per species for the complete single-cell RNA-seq preprocessing workflow, including:

- preprocessing
- quality control
- filtering
- normalization
- clustering
- doublet detection using DoubletFinder (for libraries containing more than 20,000 cells)

These scripts generate one processed Seurat object per species:

```text
OL.rds
```

The resulting files are available at GEO under accession number:

```text
GEO accession: GSE333138
```

The same scripts also contain the workflow used to isolate medulla neuroblasts, generating one:

```text
medulla_NBs.rds
```

file per species, also available at GEO under accession number:

```text
GEO accession: GSE333138
```

---

## `chose_temporal_TFs/`

This folder contains the detailed workflow for:

- trajectory analysis
- pseudotime inference
- identification of dynamically expressed transcription factors

Dynamic transcription factor expression was assessed using:

- PseudotimeDE

The main script performs trajectory analysis and temporal TF selection across species.

Additional required input data include one file for each species located in the PseudotimeDE_res_lists file:

```text
res_list.rds
```

---

## `Figure/`

This folder contains individual R scripts used to generate each figure separately.

Each script includes all analysis and visualization steps required to reproduce the corresponding figure panels.

---

## `functions/`

This folder contains different functions that are required for the execution of different parts of the analyses, including:

- LayeredFeaturePlot.R required for the representation of multiple FeaturePlots in one Feautureplot
- heatmap.R required to generate heatmaps
- plot_selection required for selecting dynamically expressed TFs, that their expression is reaching nearly 0, in at least one interval along pseudotime, and are adequately expressed.
---

## `data_geneext/`

This folder contains additional data required for figure generation, including:

- data required for Figure S3
- improved annotations for each species after GeneExt extensions

---

# Data Availability

Processed Seurat objects and additional large files are deposited in GEO:

```text
GEO accession: GSE333138
```

---

# Software Requirements

Analyses were performed in R.

Main packages include:

- Seurat
- DoubletFinder
- Slingshot
- PseudotimeDE
- tidyverse
- ggplot2

Additional package requirements are specified within individual scripts.

---

# Notes

- Each figure can be reproduced independently using the corresponding script in the `Figure/` directory.
- Large processed objects and supplementary datasets are hosted externally through GEO.
- Some analyses require intermediate files generated in previous workflow steps.
# EvoTempoPattern
