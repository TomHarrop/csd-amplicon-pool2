#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(GenomicFeatures)
library(VariantAnnotation)

txdb_file <- snakemake@input[["txdb"]]
fa_file <- snakemake@input[["fa"]]
vcf_file <- snakemake@input[["vcf"]]
tbi_file <- snakemake@input[["tbi"]]

csd_vcf <- snakemake@output[["csd"]]

# dev 
# txdb_file <- "output/005_ref/txdb.sqlite"
# fa_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
# vcf_file <- "output/035_medaka/flongle/BB24_41/round_1_phased.vcf.gz"
# tbi_file <- "output/035_medaka/flongle/BB24_41/round_1_phased.vcf.gz.tbi"
# csd_vcf <- "test/flongle-BB24_41-csd.vcf"


# read the txdb
txdb <- AnnotationDbi::loadDb(txdb_file)

# load the fasta
fa <- FaFile(fa_file, paste0(fa_file, ".fai"))

# read the vcf
tabix_file <- TabixFile(vcf_file, tbi_file)
vcf <- readVcf(tabix_file)

# get coding regions
coding <- predictCoding(vcf, txdb, fa)

# clean up / release mem
rm(vcf, fa, txdb)
gc(TRUE)

# subset by gene-level coding variants
rng <- coding[coding$GENEID == "Csd" &
                  coding$QUAL > 20 &
                  coding$CONSEQUENCE == "nonsynonymous"]
vcf_rng <- readVcf(tabix_file, param = rng)
qual(vcf_rng) <- as.integer(round(qual(vcf_rng), 0))

# write output
writeVcf(vcf_rng, csd_vcf, index = FALSE)

# log
sessionInfo()
