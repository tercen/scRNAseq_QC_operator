# scRNA-seq QC operator

##### Description
`scRNA-seq QC` performs quality-control on single-cell RNA-seq count data and returns the counts of the filtered cells and genes.

##### Usage

Input projection|.
---|---
`y-axis`        | numeric, count data, per cell 
`x-axis`        | character, cell ID
`row names`     | character, gene ID

Output relations|.
---|---
`y-axis`        | numeric, count data, per cell 
`x-axis`        | character, cell ID
`row names`     | character, gene ID

##### Details
The operator uses the QC worklfow described in the corresponding chapter of the ["Orchestrating Single-Cell Analysis"](https://osca.bioconductor.org/quality-control.html) book. For this it uses the _scRNAseq_ BioConductor package.

#### References
Amezquita, et. al. ["Orchestrating single-cell analysis with BioConductor"](https://www.nature.com/articles/s41592-019-0654-x), Nature Methods (2019)

##### See Also

#### Examples
