#!/usr/bin/env R

rm(list = ls())

metadata <- read.csv("../data/metadata/orig_metadata.csv", header=TRUE, sep=",")

# create initial popmap with all 206 samples
popmap <- metadata[c("sample", "site", "radseq_ID")]
write.table(popmap, "../data/metadata/popmap_206n_site_spp.csv", sep="\t", 
            col.names = FALSE, row.names=FALSE, quote=FALSE)
