#!/usr/bin/env R
# Takes the reads_table.tsv from the process_radtags log and assays performance

rm(list = ls())
graphics.off()
library(ggplot2)
library(stringr)
library(gridExtra)

reads_table <- read.csv("../data/process_radtags_outputs/reads_table.tsv",
                          header=TRUE, sep="\t")
reads_table <- head(reads_table,-1)  # remove control sequence

# popmap: sample, site and species
popmap <- read.csv("../data/metadata/popmap_206n_site_spp.csv",
                   header=FALSE, sep="\t")
names(popmap) <- c("Sample", "Field_site", "Species")
popmap$Species <- as.factor(popmap$Species)
popmap$Field_site <- as.factor(popmap$Field_site)
# str(popmap)

# combine process radtags data with metadata
reads_table$Sample <- str_split_fixed(reads_table$File, "[.]", 2)[,1]
reads_table <- merge(reads_table, popmap, by="Sample")


#### Overall statistics
total_reads <- sum(reads_table$Total)
total_retained <- sum(reads_table$Retained.Reads)
percent_retained <- total_retained / total_reads
print(percent_retained)

# 97.6% of reads were retained which is very high and good news. The small
# percent lost were not retained due to low quality.


#### Proportion of each sample in the total library

# calculate proportion of reads per sample in the total library
reads_table$proportion_total <- reads_table$Total / total_reads
summary(reads_table$proportion_total)
sd(reads_table$proportion_total)

# sort according to percent total
reads_table <- reads_table[order(reads_table$proportion_total),]

# ensures that plot will be ordered by percent
reads_table$File <- as.character(reads_table$File)
reads_table$File <- factor(reads_table$File, levels=reads_table$File)

# plot 
prop_lib <- ggplot(data = reads_table, aes(x = File, y = proportion_total)) +
    geom_point(stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    scale_y_continuous(expand = c(0,0.0001)) +
    labs(x = "Sample", y = "Proportion of Library")
prop_lib

# Sample 170 is contributing almost nothing to the total library. There are 
# no samples which are completely dominant.

ggsave(plot=prop_lib, 
       filename="../data/process_radtags_outputs/assay_figures/proportion_of_library.pdf",
       width=12)


# check for biases in species and field sites

# NA species = samples removed later due to being male/poor alignment
reads_no_NA <- na.omit(reads_table)

prop_lib_spp <- ggplot(data = reads_no_NA,
                       aes(x = Species, y = proportion_total)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Proportion of Library")
prop_lib_spp

prop_lib_site <- ggplot(data = reads_table,
                       aes(x = Field_site, y = proportion_total)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Field Site", y = "Proportion of Library")
prop_lib_site

# no obvious evidence of biases among any of the species or sites


#### Number of reads retained per sample

# calculate distribution of reads per sample
summary(reads_table$Retained.Reads)
sd(reads_table$Retained.Reads)

# sort by retained reads
reads_table <- reads_table[order(reads_table$Retained.Reads),]

# plot 
retained_num <- ggplot(data = reads_table, aes(x = File, y = Retained.Reads)) +
    geom_point(stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    scale_y_continuous(expand = c(0,100000)) +
    geom_hline(yintercept = 1000000, color = "red", lty=2) +
    labs(x = "Sample", y = "Retained Reads")
retained_num
# Variation in retained reads is expected. Sample 170 has an erroneously 
# low value of retained reads. 

ggsave(plot=retained_num, 
       filename="../data/process_radtags_outputs/assay_figures/num_retained_reads.pdf",
       width=12)

# As per Rivera-Colón and Catchen (2021) samples retaining fewer than 1 million
# reads could be considered too low (depending on the dataset properties), 
# here this includes sample 170 and 166. 
reads_table$File[which(reads_table$Total < 1000000)]
# When inspecting the total number of reads per sample, all apart from these 
# two samples exceed 1 million and so this seems to be a reasonable threshold 
# for this dataset.

# Rivera-Colón, A. G. and Catchen, J. (2021), Population genomics analysis with RAD, reprised: Stacks
# 2, Technical report, bioRxiv. Section: New Results Type: article.
# URL: https://www.biorxiv.org/content/10.1101/2021.11.02.466953v1


# check for biases in spp and sites

retained_num_spp <- ggplot(data = reads_no_NA, aes(x = Species, y = Retained.Reads)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Species", y = "Retained Reads")
retained_num_spp

retained_num_site <- ggplot(data = reads_table, aes(x = Field_site, y = Retained.Reads)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Field Site", y = "Retained Reads")
retained_num_site

# no obvious evidence of biases among any of the species or sites


#### Percentage of reads retained per sample

# calculate percent retained
reads_table$proportion_retained <- reads_table$Retained.Reads / reads_table$Total
summary(reads_table$proportion_retained)
sd(reads_table$proportion_retained)

# sort by retained reads percent
reads_table <- reads_table[order(reads_table$proportion_retained),]

# ensures plot ordered by percent
reads_table$File <- as.character(reads_table$File)
reads_table$File <- factor(reads_table$File, levels=reads_table$File)

retained_prop <- ggplot(data = reads_table, aes(x = File, y = proportion_retained)) +
    geom_point(stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +
    labs(x = "Sample", y = "Proportion of Reads Retained")
retained_prop

# Most samples retained a very high proportion of reads. Sample 170 again 
# had a much lower value than all others.

ggsave(plot=retained_prop, 
       filename="../data/process_radtags_outputs/assay_figures/proportion_retained_reads.pdf",
       width=12)


# spp/field site biases:

reads_no_NA <- na.omit(reads_table)
retained_prop_spp <- ggplot(data = reads_no_NA, aes(x = Species, y = proportion_retained)) +
    geom_boxplot() +
    theme_classic() +
    labs(x = "Species", y = "Proportion of Reads Retained")
retained_prop_spp

# without (ignoring) the erroneously low sample 170
retained_prop_site <- ggplot(data = reads_table, aes(x = Field_site, y = proportion_retained)) +
    geom_boxplot() +
    theme_classic() +
    ylim(0.92,0.99) +
    labs(x = "Field Site", y = "Proportion of Reads Retained")
retained_prop_site

# no obvious evidence of biases among any of the species or sites

## grid summary of plots checking for biases among spp and sites
# library(gridExtra)
spp_site_compare <- grid.arrange(prop_lib_spp, prop_lib_site,
             retained_num_spp, retained_num_site,
             retained_prop_spp, retained_prop_site)
ggsave(plot=spp_site_compare,
       filename="../data/process_radtags_outputs/assay_figures/assay_spp_site_biases.pdf",
       width=15, height=12)


#### Conclusions

# Sample 170 is an outlier in all the above assessments, it had a much lower 
# total read count than all other samples to begin and lower retained reads 
# than all others. Sample 166 also had a lower total and retained read count 
# than all other samples; according to the threshold of 1 million retained 
# reads this sample could also be rejected. 

# Samples 170 and 166 will be removed from further analyses.


## save full reads table 
write.table(reads_table, 
            file="../data/process_radtags_outputs/reads_table_assayed.tsv",
            sep = "\t", row.names = FALSE)


# save key data to be added to metadata file later
process_radtags_meta <- reads_table[c("Sample", "Retained.Reads",
                                      "proportion_retained", "proportion_total")]
names(process_radtags_meta) <- c("sample", "num_retained_reads", 
                                 "proportion_retained_reads",
                                 "proportion_library_reads")

write.csv(process_radtags_meta,
          "../data/process_radtags_outputs/process_radtags_meta.csv",
          row.names=FALSE)


## save popmap with samples 170 and 166 removed, left with 204 samples:
popmap_204n <- subset(popmap, Sample != "BB_P2_PG_170_1")
popmap_204n <- subset(popmap_204n, Sample != "BB_P3_PG_166_1")
# remove field site
popmap_204n$Field_site <- NULL

write.table(popmap_204n, file = "../data/metadata/popmap_204n_spp.csv", 
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
