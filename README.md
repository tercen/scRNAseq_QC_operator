# scRNA-seq QC

##### Description

`scRNA-seq QC` performs quality-control on single-cell RNA-seq count data and returns the counts of the filtered cells and genes.

##### Usage

Input projection|.
---|---
`y-axis`   | numeric, count data
`columns`  | character, cell ID or name
`rows`     | character, gene __name__

Output relations|.
---|---
`percent_mt`       | numeric, percent counts mapping to mitochondrial genes per cell
`nCount_RNA`       | numeric, total counts per cell
`nFeature_RNA`     | numeric, number of features (genes) detected per cell

##### Details

The operator uses the `Seurat` R package and the QC workflow described in the ["package website"](https://satijalab.org/seurat/).

##### References

> Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W. M., Zheng, S., Butler, A., ... & Satija, R. (2021). Integrated analysis of multimodal single-cell data. Cell, 184(13), 3573-3587.

[Link to Seurat reference](https://doi.org/10.1016/j.cell.2021.04.048)
