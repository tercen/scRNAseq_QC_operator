as_Seurat <- function(ctx, dim_names = "default") {
  spm <- ctx$as.sparse.matrix()
  
  if(dim_names == "default") {
    dimnames(spm) <- list(
      seq_len(ctx$rschema$nRows) - 1L,
      seq_len(ctx$cschema$nRows) - 1L
    )
  } else {
    dimnames(spm) <- list(
      ctx$rselect()[[1]],
      seq_len(ctx$cschema$nRows) - 1L
    )
  }
  CreateSeuratObject(spm)
}

seurat_to_tercen <- function(object) {
  spm <- GetAssayData(object)
  tr <- try(slot(spm, "j"))
  if(inherits(tr, "try-error")) {
    dff <- data.frame(
      i = spm@i + 1L,  
      j = as.integer(rep(1:spm@Dim[2], diff(spm@p))), 
      x = spm@x
    )
  } else {
    dff <- as.data.frame(summary(spm)) %>% 
      as_tibble() %>%
      mutate(i = as.integer(i), j = as.integer(j))
  }
  return(dff)
}