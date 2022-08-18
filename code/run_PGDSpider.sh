#!/bin/bash
# runs PGDSpider to convert files from Stacks' populations output (VCF) into the format ready for Bayescan

# requires: PGDSpider (linux install with anaconda: conda install -c bioconda pgdspider)

# activate conda env to use PGDSpider software
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ~/anaconda3/

# ruderatus - vcf
rud_spid_file=../data/PGDSpider/VCF_to_BAYESCAN_rud.spid 

# aligned to hort ref genome
rud_input=../data/ref_map_outputs/populations_vcf_rud/populations.snps.vcf
rud_output=../data/PGDSpider/rud_bayescan/pgds_vcf_bayescan_rud
PGDSpider2-cli \
    -inputfile $rud_input \
    -inputformat VCF \
    -outputfile $rud_output \
    -outputformat GESTE_BAYE_SCAN \
    -spid $rud_spid_file

# # aligned to terr ref genome
# rud_input_t=../data/ref_map_outputs/pops_vcf_rud_t/populations.snps.vcf
# rud_output_t=../data/PGDSpider/rud_bayescan/pgds_vcf_bayescan_rud_alignt
# PGDSpider2-cli \
#     -inputfile $rud_input_t \
#     -inputformat VCF \
#     -outputfile $rud_output_t \
#     -outputformat GESTE_BAYE_SCAN \
#     -spid $rud_spid_file
