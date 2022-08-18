# using VCFtools outputs to find min window size for kernel smoothing

rm(list=ls())


## Nucleotide diversity ----

## terr

terr_snps <- read.csv("../data/vcftools/terr_t_150000.snpden",
                      header=T, sep="\t")
mean(terr_snps$SNP_COUNT)
min(terr_snps$SNP_COUNT)


terr_snps <- read.csv("../data/vcftools/terr_t_150_onesigma.snpden",
                      header=T, sep="\t")
mean(terr_snps$SNP_COUNT)
min(terr_snps$SNP_COUNT)


terr_snps <- read.csv("../data/vcftools/terr_t_50000.snpden",
                      header=T, sep="\t")
mean(terr_snps$SNP_COUNT)
min(terr_snps$SNP_COUNT)



terr_vcf_pi <- read.csv("../data/vcftools/terr_t_50000.windowed.pi",
                        header=T, sep="\t")
sum(terr_vcf_pi$N_VARIANTS)

plot(terr_vcf_pi$BIN_START, terr_vcf_pi$PI)

chroms <- unique(terr_vcf_pi$CHROM)
plots <- list()

x <- terr_vcf_pi$PI
top <- quantile(x, 0.95)
bottom <- quantile(x, 0.05)

low <- subset(terr_vcf_pi, PI < bottom)

for (i in 1:length(unique(terr_vcf_pi$CHROM))) {
    # for (i in 1:1) {
    data <- terr_vcf_pi[which(terr_vcf_pi$CHROM == chroms[i]),]
    
    p <- ggplot(data, aes(x=BIN_START, y=PI)) +
        geom_point() +
        theme_classic() +
        geom_hline(yintercept=mean(data$PI, na.rm = TRUE), col="blue", lty=2) +
        geom_hline(yintercept=0, col="grey", lty=2) +
        geom_hline(yintercept = top, col="red") +
        geom_hline(yintercept = bottom, col="red") +
        labs(title=chroms[i], x="Basepair Position", y="PI")
    
    plots[[i]] <- p
}

NC_063269.1 <- plots[[1]]
NC_063270.1 <- plots[[2]]
NC_063271.1 <- plots[[3]]
NC_063272.1 <- plots[[4]]
NC_063273.1 <- plots[[5]]
NC_063274.1 <- plots[[6]]
NC_063275.1 <- plots[[7]]
NC_063276.1 <- plots[[8]]
NC_063277.1 <- plots[[9]]
NC_063278.1 <- plots[[10]]
NC_063279.1 <- plots[[11]]
NC_063280.1 <- plots[[12]]
NC_063281.1 <- plots[[13]]
NC_063282.1 <- plots[[14]]
NC_063283.1 <- plots[[15]]
NC_063284.1 <- plots[[16]]
NC_063285.1 <- plots[[17]]
NC_063286.1 <- plots[[18]]

all_terr_plots <- grid.arrange(NC_063269.1,
                               NC_063270.1,
                               NC_063271.1,
                               NC_063272.1,
                               NC_063273.1,
                               NC_063274.1,
                               NC_063275.1,
                               NC_063276.1,
                               NC_063277.1,
                               NC_063278.1,
                               NC_063279.1,
                               NC_063280.1,
                               NC_063281.1,
                               NC_063282.1,
                               NC_063283.1,
                               NC_063284.1,
                               NC_063285.1,
                               NC_063286.1,
                               nrow=4)






## B. ruderatus ----

rud_vcf_pi <- read.csv("../data/vcftools/rud_h_orig_params.windowed.pi",
                   header=TRUE, sep="\t")

sum(rud_vcf_pi$N_VARIANTS)
mean(rud_vcf_pi$N_VARIANTS)

hist(rud_vcf_pi$PI)

plot(rud_vcf_pi$BIN_START, rud_vcf_pi$PI)

ggplot(rud_vcf_pi, aes(x=BIN_START, y=PI, col=CHROM)) +
    geom_point() +
    theme_classic()


chroms <- unique(rud_vcf_pi$CHROM)
plots <- list()

x <- rud_vcf_pi$PI
top <- quantile(x, 0.95)
bottom <- quantile(x, 0.05)

for (i in 1:length(unique(rud_vcf_pi$CHROM))) {
    # for (i in 1:1) {
    data <- rud_vcf_pi[which(rud_vcf_pi$CHROM == chroms[i]),]
    
    p <- ggplot(data, aes(x=BIN_START, y=PI)) +
        geom_point() +
        theme_classic() +
        geom_hline(yintercept=mean(data$PI, na.rm = TRUE), col="blue", lty=2) +
        geom_hline(yintercept=0, col="grey", lty=2) +
        geom_hline(yintercept = top, col="red") +
        geom_hline(yintercept = bottom, col="red") +
        labs(title=chroms[i], x="Basepair Position", y="PI")
    
    plots[[i]] <- p
}


HG995188.1 <- plots[[1]]
HG995189.1 <- plots[[2]]
HG995190.1 <- plots[[3]]
HG995191.1 <- plots[[4]]
HG995192.1 <- plots[[5]]
HG995193.1 <- plots[[6]]
HG995194.1 <- plots[[7]]
HG995195.1 <- plots[[8]]
HG995196.1 <- plots[[9]]
HG995197.1 <- plots[[10]]
HG995198.1 <- plots[[11]]
HG995199.1 <- plots[[12]]
HG995200.1 <- plots[[13]]
HG995201.1 <- plots[[14]]
HG995202.1 <- plots[[15]]
HG995203.1 <- plots[[16]]
HG995204.1 <- plots[[17]]
HG995205.1 <- plots[[18]]
CAJOSO010000002.1 <- plots[[19]]
CAJOSO010000007.1 <- plots[[20]]
CAJOSO010000008.1 <- plots[[21]]

library(gridExtra)
all_chrom <- grid.arrange(HG995188.1,
                          HG995189.1,
                          HG995190.1,
                          HG995191.1,
                          HG995192.1,
                          HG995193.1,
                          HG995194.1,
                          HG995195.1,
                          HG995196.1,
                          HG995197.1,
                          HG995198.1,
                          HG995199.1,
                          HG995200.1,
                          HG995201.1,
                          HG995202.1,
                          HG995203.1,
                          HG995204.1,
                          HG995205.1, nrow=4)






unique(rud_vcf_pi$CHROM)
rud_vcf_pi$CHROM <- as.factor(rud_vcf_pi$CHROM)
rud_vcf_pi$CHROM <- factor(rud_vcf_pi$CHROM, levels=c("HG995188.1",
                                          "HG995189.1",
                                          "HG995190.1",
                                          "HG995191.1",
                                          "HG995192.1",
                                          "HG995193.1",
                                          "HG995194.1",
                                          "HG995195.1",
                                          "HG995196.1",
                                          "HG995197.1",
                                          "HG995198.1",
                                          "HG995199.1",
                                          "HG995200.1",
                                          "HG995201.1",
                                          "HG995202.1",
                                          "HG995203.1",
                                          "HG995204.1",
                                          "HG995205.1",
                                          "CAJOSO010000002.1",
                                          "CAJOSO010000007.1",
                                          "CAJOSO010000008.1"))



## from stacks smoothed

rud_pi <- rud_pi %>% 
    
    # Compute chromosome size
    group_by(Chr) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    # select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(rud_pi, ., by=c("Chr"="Chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chr, BP) %>%
    mutate( BPcumul=BP+tot)

axisdf = rud_pi %>% group_by(Chr) %>% summarize(center=( max(BPcumul) + min(BPcumul) ) / 2 )
############################################ edn stacks code





# expected pi (theta) 
n = 56
S = 3901
theta <- function(n, S) {
    denom = 0
    for (i in 1:(n-1)) {
        denom <- denom + (1/i)
    }
    return(S / denom)
}
theta(n, S)


## B. ruderatus ----

tajima_rud <- read.csv("../data/vcftools/first_test.Tajima.D",
                       header=TRUE, sep="\t")
tajima_rud <- na.omit(tajima_rud)

sum(tajima_rud$N_SNPS)

mean(tajima_rud$N_SNPS)

hist(tajima_rud$TajimaD, breaks=16)

plot(tajima_rud$BIN_START, tajima_rud$TajimaD)

ggplot(tajima_rud, aes(x=BIN_START, y=TajimaD, col=CHROM)) +
    geom_point() +
    theme_classic()

chroms <- unique(tajima_rud$CHROM)
plots <- list()

for (i in 1:length(unique(tajima_rud$CHROM))) {
# for (i in 1:1) {
    data <- tajima_rud[which(tajima_rud$CHROM == chroms[i]),]
    
    x <- data$TajimaD
    top <- quantile(x, 0.95)
    bottom <- quantile(x, 0.05)
    
    p <- ggplot(data, aes(x=BIN_START, y=TajimaD)) +
        geom_point() +
        theme_classic() +
        geom_hline(yintercept=mean(data$TajimaD, na.rm = TRUE), col="blue", lty=2) +
        geom_hline(yintercept=0, col="grey", lty=2) +
        geom_hline(yintercept = top, col="red") +
        geom_hline(yintercept = bottom, col="red") +
        labs(title=chroms[i], x="Basepair Position", y="Tajima's D")
    
    plots[[i]] <- p
}


HG995188.1 <- plots[[1]]
HG995189.1 <- plots[[2]]
HG995190.1 <- plots[[3]]
HG995191.1 <- plots[[4]]
HG995192.1 <- plots[[5]]
HG995193.1 <- plots[[6]]
HG995194.1 <- plots[[7]]
HG995195.1 <- plots[[8]]
HG995196.1 <- plots[[9]]
HG995197.1 <- plots[[10]]
HG995198.1 <- plots[[11]]
HG995199.1 <- plots[[12]]
HG995200.1 <- plots[[13]]
HG995201.1 <- plots[[14]]
HG995202.1 <- plots[[15]]
HG995203.1 <- plots[[16]]
HG995204.1 <- plots[[17]]
HG995205.1 <- plots[[18]]
CAJOSO010000002.1 <- plots[[19]]
CAJOSO010000007.1 <- plots[[20]]
CAJOSO010000008.1 <- plots[[21]]

library(gridExtra)
all_chrom <- grid.arrange(HG995188.1,
                          HG995189.1,
                          HG995190.1,
                          HG995191.1,
                          HG995192.1,
                          HG995193.1,
                          HG995194.1,
                          HG995195.1,
                          HG995196.1,
                          HG995197.1,
                          HG995198.1,
                          HG995199.1,
                          HG995200.1,
                          HG995201.1,
                          HG995202.1,
                          HG995203.1,
                          HG995204.1,
                          HG995205.1, nrow=4)
    


x <- tajima_rud$TajimaD
top <- quantile(x, 0.95)
bottom <- quantile(x, 0.05)

d <- density(tajima_rud$TajimaD)
plot(d)
abline(v=top, col="red", lty=2)
abline(v=bottom, col="red", lty=2)
abline(v=0, col="red", lty=2)

# top 5% values
top
inspect <- tajima_rud[tajima_rud$TajimaD > top,]
low <- tajima_rud[tajima_rud$TajimaD < 0,]




require(ape)
library(pegas)
rud <- read.vcf("../data/ref_map_outputs/populations_vcf_rud/populations.snps.vcf")
library(vcfR)
rud_R <- read.vcfR("../data/ref_map_outputs/populations_vcf_rud/populations.snps.vcf")

rud_DNAbin <- vcfR2DNAbin(
    rud_R,
    extract.indels = TRUE,
    consensus = FALSE,
    extract.haps = TRUE,
    unphased_as_NA = FALSE,
    asterisk_as_del = FALSE,
    ref.seq = NULL,
    start.pos = NULL,
    verbose = TRUE
)

# pegas - tajima
tajima.test(rud_DNAbin)

# sliing window
sw(rud_DNAbin, 1000, 1000, FUN = tajima.test, rowAverage = TRUE)

# # amova
# d <- dist.dna(rud_DNAbin)
# rud_popmap_subs <- read.csv("../data/metadata/popmap_rud_56n_clumpak_subpops.csv",
#                        sep="\t", header=F)
# rud_popmap <- read.csv("../data/metadata/popmap_rud_56n_site.csv",
#                        sep="\t", header=F)
# pops <- factor(rep(rud_popmap$V2, 2))
# clumpak_subpops <- factor(rep(rud_popmap_subs$V2, 2))
# amova(d ~ pops/clumpak_subpops, nperm=1000)
# amova(d ~ clumpak_subpops, nperm=10000)







hort_R <- read.vcfR("../data/ref_map_outputs/populations_vcf_hort/populations.snps.vcf")

hort_DNAbin <- vcfR2DNAbin(
    hort_R,
    extract.indels = TRUE,
    consensus = FALSE,
    extract.haps = TRUE,
    unphased_as_NA = FALSE,
    asterisk_as_del = FALSE,
    ref.seq = NULL,
    start.pos = NULL,
    verbose = TRUE
)

# pegas
tajima.test(hort_DNAbin)

# amova
d_hort <- dist.dna(hort_DNAbin)
hort_popmap <- read.csv("../data/metadata/popmap_hort_67n_site_no22.csv",
                       sep="\t", header=F)
pops_hort <- factor(rep(hort_popmap$V2, 2))
amova(d_hort ~ pops_hort, nperm=1000)




terr_R <- read.vcfR("../data/ref_map_outputs/populations_vcf_terr/populations.snps.vcf")

terr_DNAbin <- vcfR2DNAbin(
    terr_R,
    extract.indels = TRUE,
    consensus = FALSE,
    extract.haps = TRUE,
    unphased_as_NA = FALSE,
    asterisk_as_del = FALSE,
    ref.seq = NULL,
    start.pos = NULL,
    verbose = TRUE
)

# pegas
tajima.test(terr_DNAbin)




