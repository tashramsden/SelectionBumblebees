#!/usr/bin/env R

rm(list=ls())

library(coda)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(viridis)


# load the plot_bayescan function provided with Bayescan 
source("../data/bayescan/plot_R.r")

## B. RUDERATUS aligned to hortorum ref genome ----

rud1 <- read.table("../data/bayescan/rud_sites/rud_sites1_fst.txt")
rud2 <- read.table("../data/bayescan/rud_sites/rud_sites2_fst.txt")
rud3 <- read.table("../data/bayescan/rud_sites/rud_sites3_fst.txt")
rud4 <- read.table("../data/bayescan/rud_sites/rud_sites4_fst.txt")

plot_bayescan(rud1, FDR=0.05)
plot_bayescan(rud2, FDR=0.05)
plot_bayescan(rud3, FDR=0.05)
plot_bayescan(rud4, FDR=0.05)

outliers_rud1 <- rud1[which(rud1$qval < 0.05),]
outliers_rud1 <- outliers_rud1[order(outliers_rud1$qval),]
outliers_rud1

sel_files <- c("../data/bayescan/rud_sites/rud_sites1.sel",
               "../data/bayescan/rud_sites/rud_sites2.sel",
               "../data/bayescan/rud_sites/rud_sites3.sel",
               "../data/bayescan/rud_sites/rud_sites4.sel")
rud_chains <- list()

for (i in 1:4) {
    # evaluate convergence
    rud_sel <- read.table(sel_files[i], header=TRUE)
    rud_chain <- mcmc(rud_sel, thin=10)
    
    rud_chains[[i]] <- rud_chain
    
    # check the trace
    plot(rud_chain)  # hairy caterpillars
    
    # summary stats
    summary(rud_chain)
    
    # check autocorrelation
    autocorr.diag(rud_chain)
    acf(rud_sel$logL)
    
    effectiveSize(rud_chain)
}

# combined chains
combined <- mcmc.list(rud_chains)
# combined = mcmc.list(rud1_chain, rud2_chain, rud3_chain, rud4_chain)
plot(combined)

gelman.diag(combined) 
# compares within and between chain variances
# outputs = if factor = 1 -> between and within vars same = convergence
gelman.plot(combined,ask)

summary(combined)


## Higher thinning interval - more iters ----

rud1_thin50 <- read.table("../data/bayescan/rud_sites/rud_sites_thin50_1_fst.txt")
rud2_thin50 <- read.table("../data/bayescan/rud_sites/rud_sites_thin50_2_fst.txt")
rud3_thin50 <- read.table("../data/bayescan/rud_sites/rud_sites_thin50_3_fst.txt")
rud4_thin50 <- read.table("../data/bayescan/rud_sites/rud_sites_thin50_4_fst.txt")

plot_bayescan(rud1_thin50, FDR=0.05)
plot_bayescan(rud2_thin50, FDR=0.05)
plot_bayescan(rud3_thin50, FDR=0.05)
plot_bayescan(rud4_thin50, FDR=0.05)

outliers_rud1_thin50 <- rud1_thin50[which(rud1_thin50$qval < 0.01),]
outliers_rud1_thin50 <- outliers_rud1_thin50[order(outliers_rud1_thin50$qval),]
outliers_rud1_thin50

sel_files_thin50 <- c("../data/bayescan/rud_sites/rud_sites_thin50_1.sel",
                      "../data/bayescan/rud_sites/rud_sites_thin50_2.sel",
                      "../data/bayescan/rud_sites/rud_sites_thin50_3.sel",
                      "../data/bayescan/rud_sites/rud_sites_thin50_4.sel")
rud_chains_thin50 <- list()

for (i in 1:4) {
    # evaluate convergence
    rud_sel <- read.table(sel_files_thin50[i], header=TRUE)
    rud_chain <- mcmc(rud_sel, thin=50)
    
    rud_chains_thin50[[i]] <- rud_chain
    
    # check the trace
    plot(rud_chain)  # hairy caterpillars
    
    # summary stats
    print(summary(rud_chain))
    
    # check autocorrelation
    print(autocorr.diag(rud_chain))
    acf(rud_sel$Fst1)
    
    print(effectiveSize(rud_chain))
}

# combined chains
combined_thin50 <- mcmc.list(rud_chains_thin50)
plot(combined_thin50)

gelman.diag(combined_thin50) 
# compares within and between chain variances
# outputs = if factor = 1 -> between and within vars same = convergence
gelman.plot(combined_thin50, ask)

# summary stats
summary(combined_thin50)

effectiveSize(combined_thin50)


## Pr_odds 100 ----

rud1_prodds100 <- read.table("../data/bayescan/rud_sites/rud_sites_prodds100_1_fst.txt")
rud2_prodds100 <- read.table("../data/bayescan/rud_sites/rud_sites_prodds100_2_fst.txt")
rud3_prodds100 <- read.table("../data/bayescan/rud_sites/rud_sites_prodds100_3_fst.txt")
rud4_prodds100 <- read.table("../data/bayescan/rud_sites/rud_sites_prodds100_4_fst.txt")

plot_bayescan(rud1_prodds100, FDR=0.05)
plot_bayescan(rud2_prodds100, FDR=0.05)
plot_bayescan(rud3_prodds100, FDR=0.05)
plot_bayescan(rud4_prodds100, FDR=0.05)

sel_files_prodds100 <- c("../data/bayescan/rud_sites/rud_sites_prodds100_1.sel",
                         "../data/bayescan/rud_sites/rud_sites_prodds100_2.sel",
                         "../data/bayescan/rud_sites/rud_sites_prodds100_3.sel",
                         "../data/bayescan/rud_sites/rud_sites_prodds100_4.sel")
rud_chains_prodds100 <- list()

for (i in 1:4) {
    # evaluate convergence
    rud_sel <- read.table(sel_files_prodds100[i], header=TRUE)
    rud_chain <- mcmc(rud_sel, thin=10)
    
    rud_chains_prodds100[[i]] <- rud_chain
    
    # check the trace
    plot(rud_chain)  # hairy caterpillars
    
    # summary stats
    print(summary(rud_chain))
    
    # check autocorrelation
    print(autocorr.diag(rud_chain))
    acf(rud_sel$logL)
    
    print(effectiveSize(rud_chain))
}

# combined chains
combined_prodds100 <- mcmc.list(rud_chains_prodds100)
plot(combined_prodds100)

gelman.diag(combined_prodds100) 
# compares within and between chain variances
# outputs = if factor = 1 -> between and within vars same = convergence
gelman.plot(combined_prodds100, ask)

# summary stats
summary(combined_prodds100)

effectiveSize(combined_prodds100)



## OUTLIER LOCI ----

snp_ids <- read.table("../data/ref_map_outputs/populations_vcf_rud/rud_snp_ids.txt", 
                      header=F)
names(snp_ids) <- "snp_id"

# add locus ids 
rud1_thin50 <- cbind(snp_ids, rud1_thin50)

# identify outlier loci
rud1_thin50 <- rud1_thin50 %>% 
    mutate(outlier = case_when(
        qval <= 0.05 & alpha >= 0 ~ "diversifying",
        qval > 0.05 ~ "neutral",
        qval <= 0.05 & alpha < 0 ~ "balancing"
    ))
rud1_thin50$outlier <- factor(rud1_thin50$outlier)
xtabs(data=rud1_thin50, ~outlier)
# str(rud1_thin50)

locus_id <- sapply(strsplit(rud1_thin50$snp_id, ":"), "[", 1)
locus_id <- as.numeric(locus_id)
rud1_thin50 <- cbind(locus_id, rud1_thin50)

rud_outliers <- rud1_thin50[rud1_thin50$outlier=="diversifying",]
rud_outliers <- rud_outliers[order(rud_outliers$qval, decreasing = F),]
locus_ids <- c(rud_outliers$locus_id)
# str(rud_outliers)

names(rud_outliers) <- c("locus_id", "snp_id", "post_prob", "log10_PO", "qval", 
                         "alpha", "bayescan_fst", "outlier")


# get chromosome and basepair position for each outlier

rud_outliers$chromosome <- rep(NA, nrow(rud_outliers))
rud_outliers$basepair <- rep(NA, nrow(rud_outliers))

# file containing snp location info
rud_loc_file <- "../data/ref_map_outputs/populations_vcf_rud/populations.snps.vcf"

for (i in 1:nrow(rud_outliers)) {
    locus <- rud_outliers$snp_id[i]
    command <- paste("bash get_chr_bp_from_outlier.sh", locus, rud_loc_file, sep=" ")
    loc <- system(command, intern=TRUE)
    rud_outliers$chromosome[i] <- loc[1]
    rud_outliers$basepair[i] <- loc[2]
}

# save outlier info
write.table(rud_outliers, "../data/bayescan/rud_sites/outliers.csv", 
            row.names=F, quote=F)


# rud_outliers <- read.table("../data/bayescan/rud_sites/outliers.csv",
#            header=T)

# plot ruderatus bayescan results - outliers
range(rud1_thin50$qval)
# log 10 qvals 
rud1_thin50$log10_q <- log10(rud1_thin50$qval)

FDR = 0.05
signif_cutoff <- log10(FDR)

rud1_thin50_plot <- ggplot(rud1_thin50, aes(log10_q, fst)) +
    geom_point(aes(fill=outlier, color=outlier), size=2) +
    geom_vline(xintercept=signif_cutoff, linetype="dashed") +
    scale_color_manual(values=c("red", "black")) +
    labs(x = expression(paste(Log[10], "(", italic("q"), " value)")),
         y = expression(paste(italic("F")[ST]))) +
    scale_x_reverse() +
    theme_classic() +
    theme(legend.position = "none") 
rud1_thin50_plot

ggsave(plot=rud1_thin50_plot, 
       file="../data/bayescan/rud_outliers.pdf",
       width=10)
# vertical line = log10 of FDR=0.05


# exploration of alleles at outlier loci ----

rud_haplotypes <- read.csv("../data/ref_map_outputs/populations_vcf_rud/populations.haplotypes.tsv",
                           sep="\t", header=T)
names(rud_haplotypes)[names(rud_haplotypes) == "X..Catalog.Locus.ID"] <- "locus_id"

# just outliers
rud_haplotypes_outliers <- rud_haplotypes[which(rud_haplotypes$locus_id %in% locus_ids),]

# library(tidyr)
rud_haps_outliers <- gather(rud_haplotypes_outliers, 
                            sample, haplotype, 3:58)

# join geo location info from metadata
metadata <- read.csv("../data/metadata/metadata.csv", header=TRUE, sep=",")
sample_loc <- metadata[,c("sample", "site")]

rud_haps_outliers <- merge(rud_haps_outliers, sample_loc, by="sample", all.y=F, all.x=F)

# table long
rud_haps_outliers <- rud_haps_outliers %>%
    separate(haplotype, c("allele1", "allele2"), "/")
rud_haps_outliers <- gather(rud_haps_outliers, allele, alleles, c("allele1", "allele2"))

rud_haps_outliers$site <- as.factor(rud_haps_outliers$site)
rud_haps_outliers$site <- factor(rud_haps_outliers$site, levels=c("Ut", "Hd", "Os"),
                    labels=c("Upton", "Hillesden", "Ouse"))

rud_haps_outliers$alleles[rud_haps_outliers$alleles == "-"] <- NA
rud_haps_outliers <- na.omit(rud_haps_outliers)

rud_haps_outliers$locus_id <- as.factor(rud_haps_outliers$locus_id)
rud_haps_outliers$locus_id <- factor(rud_haps_outliers$locus_id, levels=c("64429",
                                                                          "1702",
                                                                          "64341",
                                                                          "1709",
                                                                          "62734",
                                                                          "22012",
                                                                          "5335"))
# save data for making heatmap in make_heatmaps.R
write.table(rud_haps_outliers,
            file="../data/bayescan/rud_sites/outliers_alleles_sites.csv",
            row.names=F, quote=F)
