#!/bin/bash
# Uses BWA to index the reference genomes
# Requires: bwa software package (http://bio-bwa.sourceforge.net/index.shtml)

# hortorum ref genome
bwa index -p ../data/align_seqs_bwa/indexed_hortorum/hortorum_index ../data/ncbi_hortorum_genome/GCA_905332935.1_iyBomHort1.1_genomic.fna.gz 
# -p = prefix for output database (will be searchable & accessed by the aligner in next step: bwa mem)
# and pass in reference genome

# Bter1.0  = old reference genome (newly annotated version 1.2 below released June so analyses re-run)
# bwa index -p ../data/align_seqs_bwa/indexed_terrestris/terrestris_index ../data/ncbi_terrestris_genome/GCA_000214255.1_Bter_1.0_genomic.fna.gz

# Bter1.2
bwa index -p ../data/align_seqs_bwa/indexed_terrestris/terrestris_index ../data/ncbi_terrestris_genome/GCF_910591885.1_iyBomTerr1.2_genomic.fna.gz
