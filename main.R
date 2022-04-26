suppressPackageStartupMessages(expr = {
  library(SingleCellExperiment)
  library(scater)
  library(tidyr)
  library(dplyr)
  library(tercen)
})

ctx = tercenCtx()

if(length(ctx$rnames) < 2) stop("At least two row factors are required.")

count_matrix <- ctx$as.matrix()

rownames(count_matrix) <- ctx$rselect()[[2]]
colnames(count_matrix) <- ctx$cselect()[[1]]

sce <- SingleCellExperiment(assays = list(counts = count_matrix))

rowData(sce)$Chr <- ctx$rselect()[[1]]

refseq <- ctx$op.value("refseq", as.character, "1")

is.mito <- which(rowData(sce)$Chr == refseq)
stats <- perCellQCMetrics(sce, subsets = list(Mito = is.mito))

high.mito <- isOutlier(stats$subsets_Mito_percent, type = "higher")
qc.libsize <- isOutlier(stats$sum, log = TRUE, type = "lower")
qc.nexprs <- isOutlier(stats$total, log = TRUE, type = "lower")
flag <- as.numeric(!(qc.libsize | qc.nexprs | high.mito))

df_out <- tibble(
  pct_counts_Mito = stats$subsets_Mito_percent,
  library_size = stats$sum,
  n_feature_detected = as.numeric(stats$detected),
  QC_flag = ifelse(flag, "pass", "fail")
) %>%
  mutate(.ci = 0:(nrow(.) - 1)) 

df_out %>%
  ctx$addNamespace() %>%
  ctx$save()

