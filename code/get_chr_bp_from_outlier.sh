#!/bin/bash
# Takes as args:
# 1: outlier locus id
# 2: path to file populations.snps.vcf 
# returns locus chromosome and basepair position

outlier=$1
file=$2

chr=$(grep $outlier $file | awk {'print $1'})
echo $chr

basepair=$(grep $outlier $file | awk {'print $2'})
echo $basepair

# chr and basepair positions can be captured in R (used in bayescan_analysis.R and run_PCAdapt.R):
# system("bash get_chr_bp_from_outlier.sh OUTLIER FILE", intern=TRUE)