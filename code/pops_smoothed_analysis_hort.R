# stacks smoothed stats + bootstrapping evaluation - B. hortorum

rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)

hort_linkage_groups <- c("HG995188.1",
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
                         "CAJOSO010000001.1",  # unplaced
                         "CAJOSO010000002.1",
                         "CAJOSO010000006.1",
                         "CAJOSO010000008.1",
                         "CAJOSO010000012.1",
                         "CAJOSO010000013.1")


## NUCLEOTIDE DIVERSITY - PI ----
# sigma 150

## h aligned ----

hort_h_sigma150 <- read.delim("../data/populations_smoothed/hort_h_smooth_sigma150/populations.sumstats.tsv",
                              skip=4)

unique(hort_h_sigma150$Chr)
hort_h_sigma150$Chr <- as.factor(hort_h_sigma150$Chr)
hort_h_sigma150$Chr <- factor(hort_h_sigma150$Chr, levels=c(hort_linkage_groups))


# remove data where smoothed pi = -1 (here overlapping loci not incl in calculation and given val -1)
hort_h_sigma150 <- subset(hort_h_sigma150, Smoothed.Pi != -1)

## remove short unplaced linkage groups
hort_h_sigma150 <- subset(hort_h_sigma150,  Chr != "CAJOSO010000001.1" &
                              Chr != "CAJOSO010000002.1" &
                              Chr != "CAJOSO010000006.1" &
                              Chr != "CAJOSO010000008.1" &
                              Chr != "CAJOSO010000012.1" &
                              Chr != "CAJOSO010000013.1")
unique(hort_h_sigma150$Chr)


# get cumulative basepair positions for plotting whole genome
hort_h_sigma150 <- hort_h_sigma150 %>% 
    
    # Compute chromosome size
    group_by(Chr) %>% 
    summarise(chr_len = max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(hort_h_sigma150, ., by=c("Chr"="Chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chr, BP) %>%
    mutate(BPcumul = BP + tot)

# get names of each chromosome (linkage group) and center of their bp for plotting x axis
axisdf = hort_h_sigma150 %>% group_by(Chr) %>% summarize(center=( max(BPcumul) + min(BPcumul) ) / 2 )

outliers_0.05 <- subset(hort_h_sigma150, Smoothed.Pi.P.value < 0.05)
outliers_0.01 <- subset(hort_h_sigma150, Smoothed.Pi.P.value < 0.01)
outliers_0.005 <- subset(hort_h_sigma150, Smoothed.Pi.P.value < 0.005)
outliers_0.001 <- subset(hort_h_sigma150, Smoothed.Pi.P.value < 0.001)
outliers_0.0005 <- subset(hort_h_sigma150, Smoothed.Pi.P.value < 0.0005)
outliers_0.0001 <- subset(hort_h_sigma150, Smoothed.Pi.P.value < 0.0001)
outliers_0.00001 <- subset(hort_h_sigma150, Smoothed.Pi.P.value < 0.00001)

length(unique(outliers_0.0001$X..Locus.ID))

write.table(outliers_0.005, file="../data/populations_smoothed/hort_h_smooth_sigma150/outliers_0.005")
write.table(outliers_0.00001, file="../data/populations_smoothed/hort_h_smooth_sigma150/outliers_0.00001")

# hist(hort_h_sigma150$Smoothed.Pi.P.value)
hist(hort_h_sigma150$Smoothed.Pi)
min(outliers_0.05$Smoothed.Pi)

top <- quantile(hort_h_sigma150$Smoothed.Pi, 0.95)
bottom <- quantile(hort_h_sigma150$Smoothed.Pi, 0.05)
low <- hort_h_sigma150[which(hort_h_sigma150$Smoothed.Pi < bottom),]

# create dataframe for plotting striped background to separate linkage groups
strip_data <- hort_h_sigma150 %>% 
    dplyr::select(BPcumul, Chr) %>% 
    group_by(Chr) %>% 
    mutate(ymin = 0, ymax = 0.009,
           xmin = min(BPcumul), xmax=max(BPcumul)) %>% 
    select(-BPcumul) %>% 
    distinct()
strip_data$fill <- rep(c("a", "b"), length.out=nrow(strip_data))

strip_data <- strip_data %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax)

sigma=150000
outliers_0.00001_pm3sigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.009,
           xmin = BPcumul - sigma*3,
           xmax = BPcumul + sigma*3) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax) %>% 
    mutate(shared = case_when(Chr == "HG995190.1" & BP <= 1575994 ~ "red",  # shared all spp
                              TRUE ~ "#fca50a"))  # not shared by all spp

outliers_0.00001_pm2sigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.009,
           xmin = BPcumul - sigma*2,
           xmax = BPcumul + sigma*2) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax) %>% 
    mutate(shared = case_when(Chr == "HG995190.1" & BP <= 1575994 ~ "red",  # shared all spp
                              TRUE ~ "#fca50a"))  # not shared by all spp


outliers_0.00001_pmsigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.009,
           xmin = BPcumul - sigma,
           xmax = BPcumul + sigma) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax) %>% 
    mutate(shared = case_when(Chr == "HG995190.1" & BP <= 1575994 ~ "red",  # shared all spp
                              TRUE ~ "#fca50a"))  # not shared by all spp


hort_h_sigma150_pi_plot <- ggplot(hort_h_sigma150, aes(x=BPcumul, y=Pi)) +
    
    geom_ribbon(data=strip_data,
                aes(y=y, xmin=xmin, xmax=xmax, group=Chr, fill=fill),
                inherit.aes=FALSE) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$Chr, breaks= axisdf$center,
                        expand=c(0.01,0.01)) +
    scale_y_continuous(expand = c(0, 0)) + 
    scale_fill_manual(name=NULL, 
                      breaks=c("a", "b"),
                      values=c("white", "grey93")) +
    
    # outlier windows, pval < 0.00001
    geom_ribbon(data=outliers_0.00001_pm3sigmawindow,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="#fca50a",
                alpha=0.2) +
    geom_ribbon(data=outliers_0.00001_pm2sigmawindow,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="#fca50a") +
    geom_ribbon(data=outliers_0.00001_pmsigmawindow,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="#fca50a") +
    
    # geom_hline(yintercept = mean(hort_h_sigma150$Smoothed.Pi), col="red") +
    geom_line(aes(group=Chr, BPcumul, Smoothed.Pi), col="grey15") +
    # geom_point(aes(group=Chr, BPcumul, Smoothed.Pi), col="grey15", size=0.2) +

    annotate(geom="point", x=44847458, y=0.0087, shape=25, col="blue", fill="blue", size=5) +
    
    # Custom the theme:
    theme_classic() +
    
    labs(y=expression(pi),
         x=expression(paste("Basepair position along ", italic("B. hortorum"), " linkage groups"))) +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 90),
          # axis.title.x = element_blank(),
          axis.text=element_text(size=20),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title=element_text(size=25))
hort_h_sigma150_pi_plot


mean(hort_h_sigma150$Smoothed.Pi)
se = sd(hort_h_sigma150$Smoothed.Pi) / sqrt(length(hort_h_sigma150$Smoothed.Pi))
se

# hort_h_outliers0.005 <- read.table("../data/populations_smoothed/hort_h_smooth_sigma150/outliers_0.005")
# hort_h_outliers0.00001 <- read.table("../data/populations_smoothed/hort_h_smooth_sigma150/outliers_0.00001")

ggsave(plot=hort_h_sigma150_pi_plot,
       file="../data/populations_smoothed/hort_h_pi.pdf",
       width=16)

ggsave(plot=hort_h_sigma150_pi_plot,
       file="../thesis/figures/hort_pi.pdf",
       width=16, height=8)

write.csv(hort_h_sigma150, 
          file="../data/populations_smoothed/hort_h_sigma150.csv",
          row.names = FALSE)

## for presentation ----

hortplot <- hort_h_sigma150_pi_plot +
    # scale_x_discrete(label=seq(1, 18, 1))
    # theme(axis.text.x = element_text(c(seq(1, 18, 1))))
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
hortplot

ggsave(plot=hortplot,
       file="../presentation/hort_pi.jpg",
       width=18, height=4)




# based on genomic regions found above - blast to terr ref genome - get flat file - contains genes within region

## get genes from biobank flat files ----
# do separately for each file - so that know which genes from which region
hort_genes_flat <- data.frame(system("bash get_genes_from_flat.sh", intern=TRUE))
names(hort_genes_flat) <- "gene"
hort_genes_flat$gene <- lapply(strsplit(hort_genes_flat$gene, '"'), "[[", 2)
hort_genes_flat


## full genes list + gene names ----
hort_h_genes <- read.csv("../data/populations_smoothed/flat_genes_hort_h/smooth_pi_hort_h_genes.csv",
                        header=T)
hort_h_genes$id  <- 1:nrow(hort_h_genes)  # to keep orig order after merge

# gene name info
hort_h_gene_names <- read.csv("../data/populations_smoothed/flat_genes_hort_h/hort_h_gene_names_DAVID.tsv",
                             header=T, sep="\t")
hort_h_gene_names$Species <- NULL
names(hort_h_gene_names) <- c("gene", "name")

hort_h_genes <- merge(hort_h_genes, hort_h_gene_names, by="gene")
hort_h_genes <- hort_h_genes[order(hort_h_genes$id),]

# reorder columns
hort_h_genes <- hort_h_genes[,c("id", "chr_hort","bp_low_h", "bp_high_h", "chr_terr", "bp_low_t", "bp_high_t", "gene", "name")]

write.csv(hort_h_genes,
          file="../data/populations_smoothed/flat_genes_hort_h/smooth_pi_hort_h_genes_with_names.csv",
          row.names=F)



# # for hort aligned - w BLAST
# get_bp_range <- function(start_bp, end_bp) {
#     current_range <- end_bp - start_bp
#     print(paste("current range: ", current_range))
#     if (current_range > 300000) {
#         wide_start <- start_bp - 300000
#         wide_end <- start_bp + 300000
#         thin_str <- paste("thin:", start_bp, end_bp, sep=" ")
#         wide_str <- paste("wide:", wide_start, wide_end, sep=" ")
#         return(paste(thin_str, wide_str, sep="    "))
#     }
#     extra <- (300000 - current_range) / 2
#     new_start <- start_bp - extra
#     new_end <- end_bp + extra
#     wide_start <- new_start - 300000
#     wide_end <- new_end + 300000
#     thin_str <- paste("thin:", new_start, new_end, sep=" ")
#     wide_str <- paste("wide:", wide_start, wide_end, sep=" ")
#     return(paste(thin_str, wide_str, sep="    "))
# }
# 
# get_bp_range(6940656, 6955947)



####### visually inspect indiv chromosomes  ----

hort_h_sigma150_pi_plot

hort_h_chr17 <- subset(hort_h_sigma150, Chr == "HG995204.1")
outliers0.00001_chr17 <- subset(outliers_0.00001, Chr == "HG995204.1")
outliers0.00001_chr17 <- outliers0.00001_chr17 %>%
    mutate(ymin=0, ymax=0.009,
           xmin = BP - sigma,
           xmax = BP + sigma) %>%
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>%
    select(-ymin_ymax)

# hort_h_chr17 <- hort_h_chr17 %>%
#     mutate(BP_plot = BP / 1000000)


hort_h_chr17_plot <- ggplot(hort_h_chr17, aes(x=BP, y=Pi)) +

    # geom_vline(xintercept = 637182, col="red") +
    # geom_vline(xintercept = 806727, col="red") +
    # geom_vline(xintercept = 721954, col="red", alpha=0.5, size=3) +

    # outlier windows, pval < 0.00001
    geom_ribbon(data=outliers0.00001_chr17,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="cyan",
                alpha=0.3) +

    geom_line(aes(group=Chr, x=BP, y=Smoothed.Pi), col="grey15") +
    geom_point(aes(x=BP, y=Smoothed.Pi), size=0.5) +

    # Custom the theme:
    theme_classic() +

    # scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6),
    #                    breaks=seq(0, max(hort_h_chr17$BP), 1000000),
    #                    expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +

    labs(y=expression(pi),
         x="Position") +
    # xlim(16000000, 18000000) +
    # geom_vline(xintercept = 18964900, col="red") +
    
    theme(legend.position="none",
          axis.text.x = element_text(),
          # axis.title.x = element_blank(),
          axis.text=element_text(size=10),
          axis.title=element_text(size=16))

hort_h_chr17_plot

