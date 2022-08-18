#!/usr/bin/env R
# Assaying performance of ref_map.pl

rm(list = ls())
library(ggplot2)
library(gridExtra)

##### HORTORUM-ALIGNED ----
## Bam statistics ----

bam_stats <- read.csv("../data/ref_map_outputs/refmap_hort_178n/bam_stats.tsv",
                        header=TRUE, sep="\t")
popmap <- read.csv("../data/metadata/popmap_178n_spp.csv", 
                   header=FALSE, sep="\t")
names(popmap) <- c("sample", "species")

# combine bam stats data with metadata
bam_stats <- merge(popmap, bam_stats, by="sample")


#### Number of primary alignments retained per sample

bam_stats <- bam_stats[order(bam_stats$primary_kept),]
# ensures that plot will be ordered
bam_stats$sample <- as.character(bam_stats$sample)
bam_stats$sample <- factor(bam_stats$sample, levels=bam_stats$sample)

num_aligns_retained <- ggplot(data = bam_stats, aes(x = sample, y = primary_kept)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    # scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Number of Retained Alignments")
num_aligns_retained

summary(bam_stats$primary_kept)

num_aligns_retained_spp <- ggplot(data = bam_stats,
       aes(x = species, y = primary_kept)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Number of Retained Alignments")
num_aligns_retained_spp

num_aligns_retained_both <- grid.arrange(num_aligns_retained, num_aligns_retained_spp, 
                                         nrow=1, widths = c(2, 1))

ggsave(plot=num_aligns_retained_both,
       filename="../data/ref_map_outputs/assay_figures/num_aligns_retained.pdf",
       width=12)

# variation in num aligns retained (was also variation in num reads per sample) - expected
# B. hort aligns best to its ref genome, followed by B. rud
# B. terr aligns least well (trend expected) - and terr aligned to hort not for subsequent analyses



#### Fraction of alignments retained per sample

bam_stats <- bam_stats[order(bam_stats$kept_frac),]
# ensures that plot will be ordered
bam_stats$sample <- as.character(bam_stats$sample)
bam_stats$sample <- factor(bam_stats$sample, levels=bam_stats$sample)

frac_aligns_retained <- ggplot(data = bam_stats, aes(x = sample, y = kept_frac)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    # scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Fraction of Retained Alignments")
frac_aligns_retained

frac_aligns_retained_spp <- ggplot(data = bam_stats,
       aes(x = species, y = kept_frac)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Fraction of Retained Alignments")
frac_aligns_retained_spp

frac_aligns_retained_both <- grid.arrange(frac_aligns_retained, frac_aligns_retained_spp, 
                                         nrow=1, widths = c(2, 1))

ggsave(plot=frac_aligns_retained_both,
       filename="../data/ref_map_outputs/assay_figures/frac_aligns_retained.pdf",
       width=12)

# most sample - majority aligns/reads kept
bam_stats[which(bam_stats$kept_frac < 0.4),"sample"]
# sample 22 (Bh) low fraction retained - remove

# again expecting frac aligns kept in B terr to be lowest 
# other 2 high


## Coverage depths ----

coverage <- read.csv("../data/ref_map_outputs/refmap_hort_178n/coverage.tsv",
                     header=TRUE, sep="\t")
coverage <- merge(popmap, coverage, by="sample")

# For mean_cov_ns, the coverage at each locus is weighted by the number of
# samples present at that locus (i.e. coverage at shared loci counts more).

# ave coverage
mean(coverage$mean_cov_ns)  # 127.1X

coverage <- coverage[order(coverage$mean_cov_ns),]
coverage$sample <- as.character(coverage$sample)
coverage$sample <- factor(coverage$sample, levels=coverage$sample)

coverage_depth <- ggplot(data = coverage, aes(x = sample, y = mean_cov_ns)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    # scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Mean Weighted Coverage Depth")
coverage_depth

coverage_depth_spp <- ggplot(data = coverage,
       aes(x = species, y = mean_cov_ns)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Mean Weighted Coverage Depth")
coverage_depth_spp


coverage_depth_both <- grid.arrange(coverage_depth, coverage_depth_spp, 
                                          nrow=1, widths = c(2, 1))

ggsave(plot=coverage_depth_both,
       filename="../data/ref_map_outputs/assay_figures/coverage_depth.pdf",
       width=12)


## Phasing rates ----

# check how well individual genotypes could be phased into two distinct
# haplotypes (diploid alleles) in each individual at each locus

phasing <- read.csv("../data/ref_map_outputs/refmap_hort_178n/phasing.tsv",
                     header=TRUE, sep="\t")
phasing <- merge(popmap, phasing, by="sample")


# number of genotyes called per sample
summary(phasing$n_gts)

phasing <- phasing[order(phasing$n_gts),]
phasing$sample <- as.character(phasing$sample)
phasing$sample <- factor(phasing$sample, levels=phasing$sample)

num_genotypes_called <- ggplot(data = phasing, aes(x = sample, y = n_gts)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    # scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Number of Genotypes Called")
num_genotypes_called

num_genotypes_called_spp <- ggplot(data = phasing,
       aes(x = species, y = n_gts)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Number of Genotypes Called")
num_genotypes_called_spp

num_genotypes_called_both <- grid.arrange(num_genotypes_called, num_genotypes_called_spp, 
                                         nrow=1, widths = c(2, 1))

ggsave(plot=num_genotypes_called_both,
       filename="../data/ref_map_outputs/assay_figures/phasing.pdf",
       width=12)


# misphasing rate - ie proportion of loci that could not be phased

summary(phasing$misphasing_rate)

phasing <- phasing[order(phasing$misphasing_rate, decreasing=TRUE),]
phasing$sample <- as.character(phasing$sample)
phasing$sample <- factor(phasing$sample, levels=phasing$sample)

misphasing <- ggplot(data = phasing, aes(x = sample, y = misphasing_rate)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    # scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Misphasing Rate")
misphasing

misphasing_spp <- ggplot(data = phasing,
                         aes(x = species, y = misphasing_rate)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Misphasing Rate")
misphasing_spp

phasing[which(phasing$misphasing_rate > 0.5),]
# sample 184 (Bh) has misphasing rate > 0.5 (0.66)

misphasing_both <- grid.arrange(misphasing, misphasing_spp, 
                                nrow=1, widths = c(2, 1))

ggsave(plot=misphasing_both,
       filename="../data/ref_map_outputs/assay_figures/misphasing_rate.pdf",
       width=12)




##### TERRESTRIS-ALIGNED ----
## Bam statistics ----

bam_stats <- read.csv("../data/ref_map_outputs/refmap_terr_178n/bam_stats.tsv",
                      header=TRUE, sep="\t")
popmap <- read.csv("../data/metadata/popmap_178n_spp.csv", 
                   header=FALSE, sep="\t")
names(popmap) <- c("sample", "species")

# combine bam stats data with metadata
bam_stats <- merge(popmap, bam_stats, by="sample")


#### Number of primary alignments retained per sample

bam_stats <- bam_stats[order(bam_stats$primary_kept),]
# ensures that plot will be ordered
bam_stats$sample <- as.character(bam_stats$sample)
bam_stats$sample <- factor(bam_stats$sample, levels=bam_stats$sample)

num_aligns_retained <- ggplot(data = bam_stats, aes(x = sample, y = primary_kept)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    # scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Number of Retained Alignments")
num_aligns_retained

summary(bam_stats$primary_kept)

num_aligns_retained_spp <- ggplot(data = bam_stats,
                                  aes(x = species, y = primary_kept)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Number of Retained Alignments")
num_aligns_retained_spp

num_aligns_retained_both <- grid.arrange(num_aligns_retained, num_aligns_retained_spp, 
                                         nrow=1, widths = c(2, 1))

ggsave(plot=num_aligns_retained_both,
       filename="../data/ref_map_outputs/assay_figures/num_aligns_retained_t.pdf",
       width=12)

#### Fraction of alignments retained per sample

bam_stats <- bam_stats[order(bam_stats$kept_frac),]
# ensures that plot will be ordered
bam_stats$sample <- as.character(bam_stats$sample)
bam_stats$sample <- factor(bam_stats$sample, levels=bam_stats$sample)

frac_aligns_retained <- ggplot(data = bam_stats, aes(x = sample, y = kept_frac)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    # scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Fraction of Retained Alignments")
frac_aligns_retained

frac_aligns_retained_spp <- ggplot(data = bam_stats,
                                   aes(x = species, y = kept_frac)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Fraction of Retained Alignments")
frac_aligns_retained_spp

frac_aligns_retained_both <- grid.arrange(frac_aligns_retained, frac_aligns_retained_spp, 
                                          nrow=1, widths = c(2, 1))

ggsave(plot=frac_aligns_retained_both,
       filename="../data/ref_map_outputs/assay_figures/frac_aligns_retained_t.pdf",
       width=12)


## Coverage depths ----

coverage <- read.csv("../data/ref_map_outputs/refmap_terr_178n/coverage.tsv",
                     header=TRUE, sep="\t")
coverage <- merge(popmap, coverage, by="sample")

# For mean_cov_ns, the coverage at each locus is weighted by the number of
# samples present at that locus (i.e. coverage at shared loci counts more).

# ave coverage
mean(coverage$mean_cov_ns)  # 128.1X

coverage <- coverage[order(coverage$mean_cov_ns),]
coverage$sample <- as.character(coverage$sample)
coverage$sample <- factor(coverage$sample, levels=coverage$sample)

coverage_depth <- ggplot(data = coverage, aes(x = sample, y = mean_cov_ns)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    # scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Mean Weighted Coverage Depth")
coverage_depth

coverage_depth_spp <- ggplot(data = coverage,
                             aes(x = species, y = mean_cov_ns)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Mean Weighted Coverage Depth")
coverage_depth_spp


coverage_depth_both <- grid.arrange(coverage_depth, coverage_depth_spp, 
                                    nrow=1, widths = c(2, 1))

ggsave(plot=coverage_depth_both,
       filename="../data/ref_map_outputs/assay_figures/coverage_depth_t.pdf",
       width=12)


## Phasing rates ----

# check how well individual genotypes could be phased into two distinct
# haplotypes (diploid alleles) in each individual at each locus

phasing <- read.csv("../data/ref_map_outputs/refmap_terr_178n/phasing.tsv",
                    header=TRUE, sep="\t")
phasing <- merge(popmap, phasing, by="sample")


# number of genotyes called per sample
summary(phasing$n_gts)

phasing <- phasing[order(phasing$n_gts),]
phasing$sample <- as.character(phasing$sample)
phasing$sample <- factor(phasing$sample, levels=phasing$sample)

num_genotypes_called <- ggplot(data = phasing, aes(x = sample, y = n_gts)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    # scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Number of Genotypes Called")
num_genotypes_called

num_genotypes_called_spp <- ggplot(data = phasing,
                                   aes(x = species, y = n_gts)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Number of Genotypes Called")
num_genotypes_called_spp

num_genotypes_called_both <- grid.arrange(num_genotypes_called, num_genotypes_called_spp, 
                                          nrow=1, widths = c(2, 1))

ggsave(plot=num_genotypes_called_both,
       filename="../data/ref_map_outputs/assay_figures/phasing_t.pdf",
       width=12)


# misphasing rate - ie proportion of loci that could not be phased

summary(phasing$misphasing_rate)

phasing <- phasing[order(phasing$misphasing_rate, decreasing=TRUE),]
phasing$sample <- as.character(phasing$sample)
phasing$sample <- factor(phasing$sample, levels=phasing$sample)

misphasing <- ggplot(data = phasing, aes(x = sample, y = misphasing_rate)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    # scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Misphasing Rate")
misphasing

misphasing_spp <- ggplot(data = phasing,
                         aes(x = species, y = misphasing_rate)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Misphasing Rate")
misphasing_spp

misphasing_both <- grid.arrange(misphasing, misphasing_spp, 
                                nrow=1, widths = c(2, 1))

ggsave(plot=misphasing_both,
       filename="../data/ref_map_outputs/assay_figures/misphasing_rate_t.pdf",
       width=12)

