library(EnsDb.Hsapiens.v75)

# fetch table with gene annotation
edb <- EnsDb.Hsapiens.v75

keys <- keys(edb, keytype="GENEID")

annots <- select(edb, keys=keys, columns=c("GENEID", "GENENAME",
                                 "SEQNAME", "SEQSTRAND",
                                 "GENESEQSTART", "GENESEQEND"),
       keytype="GENEID")

write.table(annots,
            file = "test/gene_annots.tsv",
            sep = "\t",
            quote = FALSE)
