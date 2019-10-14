#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(GenomicFeatures)

gff_file <- snakemake@input[["gff"]]
regions_file <- snakemake@output[["regions"]]

# dev
gff_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.gff"

# make a TXDB
txdb <- makeTxDbFromGFF(file = gff_file,
                        format = "gff3",
                        dataSource = "Amel_HAv3.1",
                        organism = "Apis mellifera",
                        taxonomyId = 7460)
# get features
all_cds <- cdsBy(txdb, "gene")
all_tx <- transcriptsBy(txdb, "gene")
all_exons <- exonsBy(txdb, "gene")

# which gene is that SNP in
goi <- "Csd"
goi_tx <- all_tx[[goi]]

# find the longest transcript (not run)
# goi_tx[which.max(width(goi_tx))]$tx_name

# find the longest overlapping cds
goi_cds <- all_cds[[goi]]
longest_cds <- as.data.table(goi_cds)[
    , sum(width), by = cds_name][which.max(V1), cds_name]

# extract the data for the longest cds exons
cdsoi <- goi_cds[goi_cds$cds_name == longest_cds]
cds_dt <- as.data.table(cdsoi)

# generate the samtools region format
outlines <- cds_dt[, paste0(seqnames, ":", start, "-", end)]
writeLines(outlines, regions_file, sep = " ")

# log
sessionInfo()
