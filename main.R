suppressPackageStartupMessages(expr = {
  library(Seurat)
  library(tidyr)
  library(dplyr)
  library(tercen)
})

source("./utils.R")
ctx = tercenCtx()

obj <- as_Seurat(ctx, dim_names = "rownames")
pattern <- ctx$op.value("pattern", as.character, "^OR*")
obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = pattern)

obj[[]] %>%
  select(-orig.ident) %>%
  mutate(.ci = seq_len(nrow(.)) - 1L) %>%
  ctx$addNamespace() %>%
  ctx$save()
