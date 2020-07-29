library(tercen)
library(dplyr)
library(scRNAseq)
library(SingleCellExperiment)
library(tidyr)
library(scater)
library(AnnotationDbi)
library(org.Hs.eg.db)

options("tercen.workflowId" = "7d6077b7fa4df6315a718714de00346e")
options("tercen.stepId"     = "803c9fc0-7fe0-4597-ae03-cae0ca0fa092")

getOption("tercen.workflowId")
getOption("tercen.stepId")

(ctx = tercenCtx())

count_matrix <- ctx$as.matrix()

rownames(count_matrix) <- ctx$rselect()[[1]]
colnames(count_matrix) <- ctx$cselect()[[1]]

sce <- SingleCellExperiment(assays = list(counts = count_matrix))

rowData(sce)$Symbol <- mapIds(org.Hs.eg.db,
                     keys = row.names(sce),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")


rowData(sce)$Chr <- mapIds(org.Hs.eg.db,
                              keys = row.names(sce),
                              column="CHR",
                              keytype="ENSEMBL",
                              multiVals="first")

is.mito <- which(rowData(sce)$Chr == "MT")
stats <- calculateQCMetrics(sce, subsets=list(Mito=is.mito))


