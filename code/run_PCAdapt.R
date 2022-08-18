# PCAdapt for detecting selection
# Installation: install.packages("pcadapt")

rm(list=ls())

library(pcadapt)
library(dplyr)

## B. ruderatus - aligned to hort ----

path_to_rud_h_sites_bed <- "../data/plink/rud_sites_h.bed"

# read.pcadapt to convert genotype file to pcadapt format
rud_h_sites <- read.pcadapt(path_to_rud_h_sites_bed, type="bed")


## choosing num K principal components
# default: data assumed diploid
rud_sites_h <- pcadapt(input = rud_h_sites, K=50)

# scree plot - variance explained by each PC = eigenvalues
plot(rud_sites_h, option = "screeplot")

# score plot

# sites
rud_sites_popmap <- read.csv("../data/metadata/popmap_rud_56n_site.csv",
                             sep="\t", header=F)
poplist_rud_sites <- rud_sites_popmap$V2

plot(rud_sites_h, option = "scores", pop = poplist_rud_sites)

# subpops
rud_subpops_popmap <- read.csv("../data/metadata/popmap_rud_56n_clumpak_subpops.csv",
                               sep="\t", header=F)
poplist_rud_subpops <- rud_subpops_popmap$V2

plot(rud_sites_h, option = "scores", pop = poplist_rud_subpops)
# ~ k=2 seems reasonable


## computing test stat based on the PCA

rud_sites_h <- pcadapt(rud_h_sites, K=2)
# can set param min.maf - default is 5% (pval for snps w maf < 0.05 not calculated)

summary(rud_sites_h)


## plots

# Manhattan plot
plot(rud_sites_h, option = "manhattan")  # plots log10 of pvals

# Q-Q plot
plot(rud_sites_h, option = "qqplot")  # check expected uniform distribution
# most pvals follow expected uniform dist 
# but smallest p vals smaller than expected -> outliers

# Histograms of the test-statistic and pvals
hist(rud_sites_h$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
# confirms most pvals follow ~ unif dist
# excess of small pvals = presence of outliers

plot(rud_sites_h, option = "stat.distribution")
# also see outliers


## Choosing cutoff for outlier detection

# for given alpha (0-1) snps w qval < alpha -> 
# outliers w an expected FDR bounded by alpha
# FDR = % of false +ves among list of candidate SNPs

# 3 methods, from least to most conservative:

# # q vals
# qval <- qvalue(rud_sites_h$pvalues)$qvalues
# alpha <- 0.05
# outliers <- which(qval < alpha)
# length(outliers)
# 
# Benjamini-Hochberg Procedure
# padj <- p.adjust(rud_sites_h$pvalues, method="BH")
# alpha <- 0.05
# outliers <- which(padj < alpha)
# length(outliers)
# outliers
# both of above find 83 outliers

# Bonferroni correction
padj <- p.adjust(rud_sites_h$pvalues, method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
outliers

# for k=2
# 33 outliers at alpha 0.05
# 19 at alpha 0.001
# 12 at alpha 0.0001


## association between PCs and outliers
snp_pc <- get.pc(rud_sites_h, outliers)


# check for LD bias
par(mfrow = c(2, 2))
for (i in 1:4) {
    plot(rud_sites_h$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
}
# loadings (contributions to the PCs) are evenly/randomly distributed/not clustered
# -> not determined by one genomic region (ie linkage disequilib not a problem here)


## outliers ----
## get locus IDs for outliers
rud_h_locus_ids <- read.csv("../data/ref_map_outputs/populations_vcf_rud/rud_snp_ids.txt",
                      header=F)
names(rud_h_locus_ids) <- "locus"

locus_id <- sapply(strsplit(rud_h_locus_ids$locus, ":"), "[", 1)
locus_id <- as.numeric(locus_id)
rud_h_locus_ids <- cbind(locus_id, rud_h_locus_ids)

rud_h_locus_ids <- cbind(rud_h_locus_ids, rud_sites_h$pvalues)
rud_h_locus_ids <- cbind(rud_h_locus_ids, padj)

# visually inspect outlier pvals (bonferroni adjusted)
plot(rud_h_locus_ids$locus_id, rud_h_locus_ids$padj)
abline(h=0.05, col="red")


# outliers = num (in order) of loci
rud_h_outliers <- rud_h_locus_ids[outliers,]
rud_h_outliers <- cbind(outliers, rud_h_outliers)

names(rud_h_outliers) <- c("outlier", "locus_id", "locus", "pval", "BF_padj")


# get chromosome and basepair position for each outlier

rud_h_outliers$chr <- rep(NA, length(outliers))
rud_h_outliers$basepair <- rep(NA, length(outliers))

# file containing location info
rud_h_loc_file <- "../data/ref_map_outputs/populations_vcf_rud/populations.snps.vcf"

for (i in 1:length(outliers)) {
    locus <- rud_h_outliers$locus[i]
    command <- paste("bash get_chr_bp_from_outlier.sh", locus, rud_h_loc_file, sep=" ")
    loc <- system(command, intern=TRUE)
    rud_h_outliers$chr[i] <- loc[1]
    rud_h_outliers$basepair[i] <- loc[2]
}

rud_h_outliers <- rud_h_outliers[order(rud_h_outliers$BF_padj, decreasing=T),]
rud_h_outliers$locus_id <- as.character(rud_h_outliers$locus_id)
rud_h_outliers$locus_id <- factor(rud_h_outliers$locus_id, levels=rud_h_outliers$locus_id)

ggplot(rud_h_outliers, aes(locus_id, BF_padj)) +
    geom_point() +
    theme_classic()


plot(rud_h_outliers$locus_id, rud_h_outliers$BF_padj) 

rud_h_outliers <- rud_h_outliers %>% arrange(BF_padj)

# # add extra columns for BLASTn searches (just for ease)
# rud_h_outliers$basepair <- as.numeric(rud_h_outliers$basepair)
# rud_h_outliers$low_search_BP <- rud_h_outliers$basepair - 500
# rud_h_outliers$high_search_BP <- rud_h_outliers$basepair + 500

write.table(rud_h_outliers, 
            file="../data/PCAdapt/rud_h_outliers.csv",
            row.names=F, sep=",")
