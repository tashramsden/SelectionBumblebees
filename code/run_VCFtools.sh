#!/bin/bash
# running VCFtools to calculate ave num SNPs per sliding window for different window sizes
# after running this script, use find_window_size.R to assess num SNPs based on these outputs

# test range of window sizes:
for window_size in 10000 50000 100000 150000 200000 250000;
    do

        # B. terrestris
        vcftools \
            --vcf ../data/populations_smoothed/terr_t_smooth_sigma150/populations.all.vcf \
            --SNPdensity $window_size \
            --out ../data/vcftools/terr_t_$window_size

        # B. hortorum
        vcftools \
            --vcf ../data/populations_smoothed/hort_h_smooth_sigma150/populations.all.vcf \
            --SNPdensity $window_size \
            --out ../data/vcftools/hort_h_$window_size

        # B. ruderatus
        vcftools \
            --vcf ../data/populations_smoothed/rud_h_smooth_sigma150/populations.all.vcf \
            --SNPdensity $window_size \
            --out ../data/vcftools/rud_h_$window_size

    done
