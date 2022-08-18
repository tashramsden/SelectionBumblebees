# using VCFtools outputs to find min window size for kernel smoothing
# based on num snps per window

rm(list=ls())

window_sizes = c("10000", "50000", "100000", "150000", "200000", "250000")

## terr
mean_snps_per_window_t <- c()
for (i in window_sizes) {
    file_name = paste("terr_t_", i, ".snpden", sep="")
    snp_file = paste("../data/vcftools/", file_name, sep="")
    terr_snps <- read.csv(snp_file, header=T, sep="\t")
    mean_snps <- mean(terr_snps$SNP_COUNT)
    mean_snps_per_window_t <- c(mean_snps_per_window_t, terr_snps)
    print(mean_snps)
}

## hort
mean_snps_per_window_h <- c()
for (i in window_sizes) {
    file_name = paste("hort_h_", i, ".snpden", sep="")
    snp_file = paste("../data/vcftools/", file_name, sep="")
    terr_snps <- read.csv(snp_file, header=T, sep="\t")
    mean_snps <- mean(terr_snps$SNP_COUNT)
    mean_snps_per_window_h <- c(mean_snps_per_window_h, terr_snps)
    print(mean_snps)
}

## rud
mean_snps_per_window_r <- c()
for (i in window_sizes) {
    file_name = paste("rud_h_", i, ".snpden", sep="")
    snp_file = paste("../data/vcftools/", file_name, sep="")
    terr_snps <- read.csv(snp_file, header=T, sep="\t")
    mean_snps <- mean(terr_snps$SNP_COUNT)
    mean_snps_per_window_r <- c(mean_snps_per_window_r, terr_snps)
    print(mean_snps)
}

## for kernel smoothing use minimum sigma 100000 to ensure enough SNPs per window (based on pers comms w Stacks authors)

