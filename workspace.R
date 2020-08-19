library(scRNAseq)
library(SingleCellExperiment)
library(scater)
library(tidyr)
library(dplyr)
library(tercen)

options("tercen.workflowId" = "7d6077b7fa4df6315a718714de00346e")
options("tercen.stepId"     = "ceff76f8-85b3-4738-a873-4d5fad13b345")

getOption("tercen.workflowId")
getOption("tercen.stepId")

(ctx = tercenCtx())

count_matrix <- ctx$as.matrix()

rownames(count_matrix) <- ctx$rselect()[[2]]
colnames(count_matrix) <- ctx$cselect()[[1]]

sce <- SingleCellExperiment(assays = list(counts = count_matrix))

rowData(sce)$Chr <- ctx$rselect()[[1]]

is.mito <- which(rowData(sce)$Chr == "MT")
stats <- calculateQCMetrics(sce, feature_controls = list(Mito=is.mito))
high.mito <- isOutlier(stats$pct_counts_Mito, type="higher")

qc.libsize <- isOutlier(stats$total_counts, log=TRUE, type="lower")
qc.nexprs <- isOutlier(stats$total_features_by_counts, log=TRUE, type="lower")

output <- tibble(pct_counts_Mito = stats$pct_counts_Mito,
                 library_size = stats$total_counts,
                 n_feature_detected = stats$total_features_by_counts,
                 passes_QC = !(qc.libsize | qc.nexprs | high.mito))
output$.ci <- 0:(nrow(output)-1)

row_indexes <- ctx %>% select(.ri, .ci)
output <- left_join(row_indexes, output)

ctx$addNamespace(output) %>%
  ctx$save()
