#!/usr/bin/env R
## Before running ref_map collate all metadata and add to it which samples will be removed
# create final popmap w/o removed samples

rm(list = ls())
library(stringr)

metadata <- read.csv("../data/metadata/orig_metadata.csv", header=TRUE, sep=",")


# add process radtags results metadata
process_radtags_meta <- read.csv("../data/process_radtags_outputs/process_radtags_meta.csv", 
                                 header = TRUE)
metadata_new <- merge(metadata, process_radtags_meta, by="sample")

# add bwa alignment metadata
hort_align_meta <- read.csv("../data/align_seqs_bwa/hort_align_meta.csv",
                            header = TRUE)
terr_align_meta <- read.csv("../data/align_seqs_bwa/terr_align_meta.csv",
                            header = TRUE)

metadata_new <- merge(metadata_new, hort_align_meta, by="sample", all=TRUE)
metadata_new <- merge(metadata_new, terr_align_meta, by="sample", all=TRUE)



# siblings ----
# identify sibling groups and only keep one per group w highest num retained reads

siblings <- metadata_new[which(!is.na(metadata_new$sibling_sets)),]
num_siblings <- nrow(siblings)
# 29 siblings

siblings$family <- gsub("^.*?_", "", siblings$sibling_sets)
siblings$family <- as.factor(siblings$family)

families <- unique(siblings$family)
num_families <- length(families)
# 12 different families

# only keep one indiv from each family - w highest num retained reads
sibs_to_keep <- rep(NA, num_families)
all_sibs_to_remove <- c()  # get rid of all siblings other than one from each fam

for (i in 1:num_families) {
    family <- siblings[siblings$family == families[i],]
    sib_to_keep <- family[which.max(family$num_retained_reads),]
    sibs_to_keep[i] <- sib_to_keep$index
    sibs_to_remove <- family[family$index != sib_to_keep$index,]
    all_sibs_to_remove <- c(all_sibs_to_remove, sibs_to_remove$index)
}

print(all_sibs_to_remove)
length(all_sibs_to_remove)
# 17 siblings will be removed

print(sibs_to_keep)
# 1 from each family w highest retained reads kept (12 kept)

all_sibs <- c(sibs_to_keep, all_sibs_to_remove)
sib_df <- data.frame(all_sibs)
sib_df$keep <- rep(NA, length(all_sibs))

for (i in 1:length(all_sibs)) {
    if (all_sibs[i] %in% sibs_to_keep) {
        sib_df$keep[i] <- TRUE
    } else {
        sib_df$keep[i] <- FALSE
    }
}

names(sib_df) <- c("index", "keep_sibling") 
fam_groups <- siblings[c("index", "family")]
sibling_info <- merge(fam_groups, sib_df, by="index")

## add sibling info to metadata
metadata_new <- merge(metadata_new, sibling_info, by="index", all=TRUE)


# samples to be removed from analyses ----

samples_to_remove <- c()

# low retained reads
samples_to_remove <- c(samples_to_remove, 170, 166, 158)

# poor alignment to reference genomes
samples_to_remove <- c(samples_to_remove, 205, 237, 242, 22)

# males
samples_to_remove <- c(samples_to_remove, 24, 64, 68)

# low F value
samples_to_remove <- c(samples_to_remove, 124)

# siblings with not highest num retained reads of family
samples_to_remove <- c(samples_to_remove, all_sibs_to_remove)

length(samples_to_remove)

index <- c(metadata_new$index)
remove_df <- data.frame(index) 
remove_df$remove_sample <- rep(NA, 206)

for (i in 1:length(index)) {
    if (index[i] %in% samples_to_remove) {
        remove_df$remove_sample[i] <- TRUE
    } else {
        remove_df$remove_sample[i] <- FALSE
    }
}

sum(remove_df$remove_sample==TRUE)
sum(remove_df$remove_sample==FALSE)
# keeping 178 samples, removing 28

metadata_new <- merge(metadata_new, remove_df, by="index", all=TRUE)

# also prep clumpak subpops categories
metadata_new$cp_subpop <- str_split_fixed(metadata_new$clumpak_subpop, "_", 3)[,2]

write.csv(metadata_new, "../data/metadata/metadata.csv", row.names = FALSE)



## create final popmap w 178 samples

metadata_178 <- metadata_new[!metadata_new$remove_sample == TRUE,]
table(metadata_178$radseq_ID)
# 67 Bh, 56 Br, 55 Bt

popmap_final <- metadata_178[c("sample", "radseq_ID")]

write.table(popmap_final, "../data/metadata/popmap_178n_spp.csv", sep="\t", 
            col.names = FALSE, row.names=FALSE, quote=FALSE)


# B. hortorum
metadata_hort <- metadata_178[metadata_178$radseq_ID == "Bh",]
popmap_hort_sites <- metadata_hort[c("sample", "site")]
write.table(popmap_hort_sites, "../data/metadata/popmap_hort_67n_site.csv", sep="\t", 
            col.names = FALSE, row.names=FALSE, quote=FALSE)
table(popmap_hort_sites$site)


# B. ruderatus
metadata_rud <- metadata_178[metadata_178$radseq_ID == "Br",]
popmap_rud_sites <- metadata_rud[c("sample", "site")]
write.table(popmap_rud_sites, "../data/metadata/popmap_rud_56n_site.csv", sep="\t", 
            col.names = FALSE, row.names=FALSE, quote=FALSE)
table(popmap_rud_sites$site)


# B. terrestris
metadata_terr <- metadata_178[metadata_178$radseq_ID == "Bt",]
popmap_terr_sites <- metadata_terr[c("sample", "site")]
write.table(popmap_terr_sites, "../data/metadata/popmap_terr_55n_site.csv", sep="\t",
            col.names = FALSE, row.names=FALSE, quote=FALSE)
table(popmap_terr_sites$site)
