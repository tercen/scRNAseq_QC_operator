library(scRNAseq)
library(SingleCellExperiment)
library(scater)
library(tidyr)
library(dplyr)
library(tercen)

options("tercen.workflowId" = "4be77e8b64b21d32888498101300ee68")
options("tercen.stepId"     = "13752223-2f8a-4e02-839a-884d386bbb56")

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
                 n_feature_detected = as.numeric(stats$total_features_by_counts))#,
            #     passes_QC = !(qc.libsize | qc.nexprs | high.mito))
output$.ci <- (0:(nrow(output)-1))

full_schema <- ctx %>% select(.ci, .ri) %>%
  left_join(output)

#output <- bind_cols(ctx %>% select(.ci, .ri) %>% dplyr::filter(.ri == 0) %>% select(.ci),
#                    output)

#output$.ci <- as.double(output$.ci)

ctx$addNamespace(full_schema) %>%
  ctx$save()
