#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(AnnotationDbi)
library(GenomicFeatures)

gff_file <- snakemake@input[["gff"]]
txdb_file <- snakemake@output[["txdb"]]

# read the txdb
txdb <- makeTxDbFromGFF(file = gff_file,
                        format = "gff3",
                        dataSource = "Amel_HAv3.1",
                        organism = "Apis mellifera",
                        taxonomyId = 7460)

# write output
AnnotationDbi::saveDb(txdb, txdb_file)

# log
sessionInfo()
