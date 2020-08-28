library(scRNAseq)
library(SingleCellExperiment)
library(scater)
library(tidyr)
library(dplyr)
library(tercen)

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
                 n_feature_detected = as.numeric(stats$total_features_by_counts))#,
output$.ci <- as.double(0:(nrow(output)-1))

ctx$addNamespace(output) %>%
  ctx$save()
