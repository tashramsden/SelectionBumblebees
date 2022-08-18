#!/usr/bin/env R
# Evaluates the alignments of each sample to the reference genomes (B. terrestris and B. hortorum)

rm(list = ls())
graphics.off()
library(ggplot2)
library(dplyr)

# B. hortorum alignment stats
percent_align_hort <- read.csv("../data/align_seqs_bwa/percent_mapped_to_hortorum.csv",
                               header=FALSE, sep=" ")
names(percent_align_hort) <- c("Sample", "Percent_mapped", "Primary_perc_mapped")

# B. terrestris alignement stats
percent_align_terr <- read.csv("../data/align_seqs_bwa/percent_mapped_to_terrestris.csv",
                               header=FALSE, sep=" ")
names(percent_align_terr) <- c("Sample", "Percent_mapped", "Primary_perc_mapped")
str(percent_align_terr)

setdiff(percent_align_hort$Sample, percent_align_terr$Sample)  # all samples present

# remove % symbol from percentages 
percent_align_hort$Percent_mapped <- gsub("%$", "", percent_align_hort$Percent_mapped)
percent_align_hort$Percent_mapped <- as.numeric(percent_align_hort$Percent_mapped)
str(percent_align_hort)

percent_align_terr$Percent_mapped <- gsub("%$", "", percent_align_terr$Percent_mapped)
percent_align_terr$Percent_mapped <- as.numeric(percent_align_terr$Percent_mapped)
str(percent_align_terr)


# metadata: sample, site and species
popmap <- read.csv("../data/metadata/popmap_204n_spp.csv", 
                   header=FALSE, sep="\t")
names(popmap) <- c("Sample", "Species")
popmap$Species <- as.factor(popmap$Species)

# combine spp ID with sample alignment info
hort_align <- merge(popmap, percent_align_hort, by="Sample")
terr_align <- merge(popmap, percent_align_terr, by="Sample")


## B. hortorum alignment
hist(hort_align$Percent_mapped)

hort_align_plot <- ggplot(data = hort_align, aes(x = Species, y = Percent_mapped)) +
    geom_point(position="jitter") +
    theme_classic() +
    labs(x = "Species", 
         y = expression(paste("Alignment Percentage to ", italic("B. hortorum"), " Reference Genome"))) +
    geom_hline(yintercept = 60, color = "red", lty=2)
hort_align_plot

ggsave(plot=hort_align_plot,
       filename="../data/align_seqs_bwa/assay_figures/align_to_hort_all.pdf",
       width=12)

hort_align[which(hort_align$Percent_mapped < 60),"Sample"]

# without the samples w no spp assignment (will be removed anyway)
hort_align_no_NA <- na.omit(hort_align)
hort_align_plot2 <- ggplot(data = hort_align_no_NA, aes(x = Species, y = Percent_mapped)) +
    geom_point(position="jitter") +
    theme_classic() +
    labs(x = "Species", 
         y = expression(paste("Alignment Percentage to ", italic("B. hortorum"), " Reference Genome"))) +
    geom_hline(yintercept = 60, color = "red", lty=2)
hort_align_plot2

ggsave(plot=hort_align_plot2,
       filename="../data/align_seqs_bwa/assay_figures/align_to_hort_no_NA.pdf",
       width=12)

# w/o terr samples - these not aligned to hort ref for analyses
hort_align_only_hort_rud <- subset(hort_align_no_NA, Species == "Bh" | Species == "Br")
hort_align_plot3 <- ggplot(data = hort_align_only_hort_rud, aes(x = Species, y = Percent_mapped)) +
    geom_point(position="jitter") +
    theme_classic() +
    labs(x = "Species", 
         y = expression(paste("Alignment Percentage to ", italic("B. hortorum"), " Reference Genome"))) +
    geom_hline(yintercept = 60, color = "red", lty=2)
hort_align_plot3

ggsave(plot=hort_align_plot3,
       filename="../data/align_seqs_bwa/assay_figures/align_to_hort.pdf",
       width=12)

hort_align %>% 
    group_by(Species) %>% 
    summarize(mean = mean(Percent_mapped),
              max = max(Percent_mapped),
              min = min(Percent_mapped))


hort_align_only_hort_rud[which(hort_align_only_hort_rud$Percent_mapped < 60),]
to_remove <- hort_align_only_hort_rud[which(hort_align_only_hort_rud$Percent_mapped < 60),"Sample"]
# 3 samples which align to < 40% of ref genome: 205, 237, 242 (NA - from field 
# ID these were B. terrestris but will be removed from samples - v poor align)
# 4 < 60%, above + 22 (Bh)

## only samples kept - summary stats of alignment:
final_hort_align_samples <- subset(hort_align_only_hort_rud, Sample != to_remove)

tab <- final_hort_align_samples %>% 
    group_by(Species) %>% 
    summarize(mean = mean(Percent_mapped),
              max = max(Percent_mapped),
              min = min(Percent_mapped),
              # standard error
              SE = sd(Percent_mapped) / sqrt(length(Percent_mapped)) )
as.data.frame(tab)



## B. terrestris alignment ----
hist(terr_align$Percent_mapped)
hist(terr_align_OLD$Percent_mapped)


terr_align_plot <- ggplot(data = terr_align, aes(x = Species, y = Percent_mapped)) +
    geom_point(position="jitter") +
    theme_classic() +
    labs(x = "Species", 
         y = expression(paste("Alignment Percentage to ", italic("B. terrestris"), " Reference Genome"))) +
    ylim(19, 100) +
    geom_hline(yintercept = 60, color = "red", lty = 2)
terr_align_plot

ggsave(plot=terr_align_plot,
       filename="../data/align_seqs_bwa/assay_figures/align_to_terr_all.pdf",
       width=12)

# without the samples w no spp assignment (will be removed anyway)
terr_align_no_NA <- na.omit(terr_align)
terr_align_plot2 <- ggplot(data = terr_align_no_NA, aes(x = Species, y = Percent_mapped)) +
    geom_point(position="jitter") +
    theme_classic() +
    labs(x = "Species", 
         y = expression(paste("Alignment Percentage to ", italic("B. terrestris"), " Reference Genome"))) +
    ylim(19, 100) +
    geom_hline(yintercept = 50, color = "red", lty = 2)
terr_align_plot2

ggsave(plot=terr_align_plot2,
       filename="../data/align_seqs_bwa/assay_figures/align_to_terr_no_NA.pdf",
       width=12)

terr_align_plot2 <- ggplot(data = terr_align_no_NA, aes(x = Species, y = Percent_mapped)) +
    geom_point(position="jitter") +
    theme_classic() +
    labs(x = "Species", 
         y = expression(paste("Alignment Percentage to ", italic("B. terrestris"), " Reference Genome"))) +
    ylim(19, 100) +
    # cut-off of 50% for hort and rud since these more distant relatives of terr so poorer alignment expected
    geom_segment(aes(x=0.5, xend=2.5, y=50, yend=50), lty=2, col="red") +
    # terr cut-off 60% - same as cutoff for hort and rud to hort ref genome
    geom_segment(aes(x=2.5, xend=3.5, y=60, yend=60), lty=2, col="red")
terr_align_plot2

ggsave(plot=terr_align_plot2,
       filename="../data/align_seqs_bwa/assay_figures/align_to_terr_diff_cutoffs.pdf",
       width=12)

terr_align %>% 
    group_by(Species) %>% 
    summarize(mean = mean(Percent_mapped),
              max = max(Percent_mapped),
              min = min(Percent_mapped))


terr_align[which(terr_align$Percent_mapped < 60),]
# 4 samples which align to < 40% of ref genome: 22 (Bh), 205, 237, 242 (NA - Bt from field ID)
# 5 < 50%, above + 162 (Br)

# remove samples aligning < 60%
terr_align_no_NA[which(terr_align_no_NA$Percent_mapped < 50),]
to_remove_terr <- terr_align_no_NA[which(terr_align_no_NA$Percent_mapped < 50),"Sample"]

## only samples kept - summary stats of alignment:
final_terr_align_samples <- subset(terr_align_no_NA, Sample != to_remove)

table(final_terr_align_samples$Species)

tab <- final_terr_align_samples %>% 
    group_by(Species) %>% 
    summarize(mean = mean(Percent_mapped),
              max = max(Percent_mapped),
              min = min(Percent_mapped),
              # standard error
              SE = sd(Percent_mapped) / sqrt(length(Percent_mapped)) )
as.data.frame(tab)

## Check stat differences between aligning each spp to each ref genome:

# B . terrestris
terr_to_terr <- final_terr_align_samples[final_terr_align_samples$Species == "Bt",]
summary(terr_to_terr$Percent_mapped)

# terr_to_hort <- final_hort_align_samples[final_hort_align_samples$Species == "Bt",]
# summary(terr_to_hort$Percent_mapped)
# t.test(terr_to_terr$Percent_mapped, terr_to_hort$Percent_mapped)
# B. terr aligned to its own ref genome is sig higher than alignment to B. hort ref genome 


# B. hortorum
hort_to_terr <- final_terr_align_samples[final_terr_align_samples$Species == "Bh",]
summary(hort_to_terr$Percent_mapped)

hort_to_hort <- final_hort_align_samples[final_hort_align_samples$Species == "Bh",]
summary(hort_to_hort$Percent_mapped)

t.test(hort_to_terr$Percent_mapped, hort_to_hort$Percent_mapped)


# B. ruderatus
rud_to_terr <- final_terr_align_samples[final_terr_align_samples$Species == "Br",]
summary(rud_to_terr$Percent_mapped)

rud_to_hort <- final_hort_align_samples[final_hort_align_samples$Species == "Br",]
summary(rud_to_hort$Percent_mapped)

t.test(rud_to_terr$Percent_mapped, rud_to_hort$Percent_mapped)

## B. hort and rud align sig better to B. hort genome than B. terr


## test if diff between hort and rud aligning to hort genome
t.test(rud_to_hort$Percent_mapped, hort_to_hort$Percent_mapped)


## Conclusions 

# samples 205, 237 and 242 do not align to either reference genome successfully
# (v low alignments) they will be removed from subsequent analyses
# sample 22 also poor alignment to hort and terr ref genomes - remove

# clear that B. hortorum and B. ruderatus samples align better to the B. hortorum reference genome
# they both align very successfully to this genome ~ 97%
# (expected - they are closely related species)

# B. terrestris aligns well to both reference genomes - more to its own


## save key stats to be added to metadata file later
meta_hort <- hort_align[c("Sample", "Percent_mapped")]
names(meta_hort) <- c("sample", "hort_align_percent")

meta_terr <- terr_align[c("Sample", "Percent_mapped")]
names(meta_terr) <- c("sample", "terr_align_percent")

write.csv(meta_hort, "../data/align_seqs_bwa/hort_align_meta.csv", row.names=FALSE)
write.csv(meta_terr, "../data/align_seqs_bwa/terr_align_meta.csv", row.names=FALSE)

