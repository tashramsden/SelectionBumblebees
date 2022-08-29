# stacks smoothed stats + bootstrapping evaluation - B. ruderatus

rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(grDevices)

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
                         "CAJOSO010000002.1",  # unplaced
                         "CAJOSO010000007.1",
                         "CAJOSO010000008.1")


## NUCLEOTIDE DIVERSITY - PI ----
# sigma 150

## h aligned ----

rud_h_sigma150 <- read.delim("../data/populations_smoothed/rud_h_smooth_sigma150/populations.sumstats.tsv",
                              skip=3)

unique(rud_h_sigma150$Chr)
rud_h_sigma150$Chr <- as.factor(rud_h_sigma150$Chr)
rud_h_sigma150$Chr <- factor(rud_h_sigma150$Chr, levels=c(hort_linkage_groups))


# remove data where smoothed pi = -1 (here overlapping loci not incl in calculation and given val -1)
rud_h_sigma150 <- subset(rud_h_sigma150, Smoothed.Pi != -1)

## remove srud unplaced linkage groups
rud_h_sigma150 <- subset(rud_h_sigma150,Chr != "CAJOSO010000002.1" &
                              Chr != "CAJOSO010000007.1" &
                              Chr != "CAJOSO010000008.1")
unique(rud_h_sigma150$Chr)


# get cumulative basepair positions for plotting whole genome
rud_h_sigma150 <- rud_h_sigma150 %>% 
    
    # Compute chromosome size
    group_by(Chr) %>% 
    summarise(chr_len = max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(rud_h_sigma150, ., by=c("Chr"="Chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chr, BP) %>%
    mutate(BPcumul = BP + tot)

# get names of each chromosome (linkage group) and center of their bp for plotting x axis
axisdf = rud_h_sigma150 %>% group_by(Chr) %>% summarize(center=( max(BPcumul) + min(BPcumul) ) / 2 )

outliers_0.05 <- subset(rud_h_sigma150, Smoothed.Pi.P.value < 0.05)
outliers_0.01 <- subset(rud_h_sigma150, Smoothed.Pi.P.value < 0.01)
outliers_0.005 <- subset(rud_h_sigma150, Smoothed.Pi.P.value < 0.005)
outliers_0.001 <- subset(rud_h_sigma150, Smoothed.Pi.P.value < 0.001)
outliers_0.0005 <- subset(rud_h_sigma150, Smoothed.Pi.P.value < 0.0005)
outliers_0.0001 <- subset(rud_h_sigma150, Smoothed.Pi.P.value < 0.0001)
outliers_0.00001 <- subset(rud_h_sigma150, Smoothed.Pi.P.value < 0.00001)

length(unique(outliers_0.0001$X..Locus.ID))

write.table(outliers_0.005, file="../data/populations_smoothed/rud_h_smooth_sigma150/outliers_0.005")
write.table(outliers_0.00001, file="../data/populations_smoothed/rud_h_smooth_sigma150/outliers_0.00001")

# hist(rud_h_sigma150$Smoothed.Pi.P.value)
hist(rud_h_sigma150$Smoothed.Pi)
min(outliers_0.05$Smoothed.Pi)

top <- quantile(rud_h_sigma150$Smoothed.Pi, 0.95)
bottom <- quantile(rud_h_sigma150$Smoothed.Pi, 0.05)
low <- rud_h_sigma150[which(rud_h_sigma150$Smoothed.Pi < bottom),]

# create dataframe for plotting striped background to separate linkage groups
strip_data <- rud_h_sigma150 %>% 
    dplyr::select(BPcumul, Chr) %>% 
    group_by(Chr) %>% 
    mutate(ymin = 0, ymax = 0.006,
           xmin = min(BPcumul), xmax=max(BPcumul)) %>% 
    select(-BPcumul) %>% 
    distinct()
strip_data$fill <- rep(c("a", "b"), length.out=nrow(strip_data))

strip_data <- strip_data %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax)

sigma=150000
outliers_0.00001_pm3sigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.006,
           xmin = BPcumul - sigma*3,
           xmax = BPcumul + sigma*3) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax) %>% 
    mutate(shared = case_when(Chr == "HG995190.1" & BP <= 1124565 ~ "red",  # shared all spp
                              TRUE ~ "#fca50a"))  # not shared by all spp

outliers_0.00001_pm2sigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.006,
           xmin = BPcumul - sigma*2,
           xmax = BPcumul + sigma*2) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax) %>% 
    mutate(shared = case_when(Chr == "HG995190.1" & BP <= 1124565 ~ "red",  # shared all spp
                              TRUE ~ "#fca50a"))  # not shared by all spp

outliers_0.00001_pmsigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.006,
           xmin = BPcumul - sigma,
           xmax = BPcumul + sigma) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax) %>% 
    mutate(shared = case_when(Chr == "HG995190.1" & BP <= 1124565 ~ "red",  # shared all spp
                              TRUE ~ "#fca50a"))  # not shared by all spp

rud_h_sigma150_pi_plot <- ggplot(rud_h_sigma150, aes(x=BPcumul, y=Pi)) +
    
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
                      # values=alpha(c("white", "#f0f921"), 0.4)) +
                      # alpha= alpha(c(1, 0.1), 0.2)) + # fcffa4
                      # values=c("white", "#DFEAF9")) +
    
    
    # outlier windows, pval < 0.00001
    geom_ribbon(data=outliers_0.00001_pm3sigmawindow,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="#fca50a") +
    geom_ribbon(data=outliers_0.00001_pm2sigmawindow,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="#fca50a") +
    geom_ribbon(data=outliers_0.00001_pmsigmawindow,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="#fca50a") +
    
    # geom_hline(yintercept = mean(rud_h_sigma150$Smoothed.Pi), col="red") +
    
    geom_line(aes(group=Chr, BPcumul, Smoothed.Pi), col="grey15") +
    # geom_point(aes(group=Chr, BPcumul, Smoothed.Pi), col="grey15", size=0.2) +
    
    # outlier markers
    # annotate(geom="text", x=c(outliers_0.005$BPcumul), y=0.0087, label="*", col="blue") +
    # annotate(geom="text", x=c(outliers_0.00001$BPcumul), y=0.009, label="*", col="red") +
    
    annotate(geom="point", x=44485844, y=0.0057, shape=25, col="blue", fill="blue", size=5) +
  
    # Custom the theme:
    theme_classic() +
    
    labs(y=expression(pi),
         x=expression(paste("Basepair position along ", italic("B. hortorum"), " linkage groups"))) +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 90),
          # axis.title.x = element_blank(),
          axis.text=element_text(size=20),
          axis.title=element_text(size=25),
          plot.margin = margin(t = 12, r = 5.5, b = 5.5, l = 5.5, unit = "pt"))

rud_h_sigma150_pi_plot

rud_h_pi_plot_for_combo <- rud_h_sigma150_pi_plot + 
    theme(legend.position="none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text=element_text(size=20),
          axis.title=element_text(size=25))
rud_h_pi_plot_for_combo


mean(rud_h_sigma150$Smoothed.Pi)
se = sd(rud_h_sigma150$Smoothed.Pi) / sqrt(length(rud_h_sigma150$Smoothed.Pi))
se

# rud_h_outliers0.005 <- read.table("../data/populations_smoothed/rud_h_smooth_sigma150/outliers_0.005")
# rud_h_outliers0.00001 <- read.table("../data/populations_smoothed/rud_h_smooth_sigma150/outliers_0.00001")

ggsave(plot=rud_h_sigma150_pi_plot,
       file="../data/populations_smoothed/rud_h_smooth_pi.pdf",
       width=16)

# ggsave(plot=rud_h_sigma150_pi_plot,
#        file="../thesis/figures/rud_pi.pdf",
#        width=16)

write.csv(rud_h_sigma150, 
          file="../data/populations_smoothed/rud_h_sigma150.csv",
          row.names = FALSE)

## for presentation ----

rudplot <- rud_h_sigma150_pi_plot +
  # scale_x_discrete(label=seq(1, 18, 1))
  # theme(axis.text.x = element_text(c(seq(1, 18, 1))))
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
rudplot

ggsave(plot=rudplot,
       file="../presentation/rud_pi.jpg",
       width=18, height=4)





## find gene bp search ----
# for hort aligned - w BLAST
get_bp_range <- function(start_bp, end_bp) {
    current_range <- end_bp - start_bp
    print(paste("current range: ", current_range))
    if (current_range > 300000) {
        wide_start <- start_bp - 300000
        wide_end <- start_bp + 300000
        thin_str <- paste("thin:", start_bp, end_bp, sep=" ")
        wide_str <- paste("wide:", wide_start, wide_end, sep=" ")
        return(paste(thin_str, wide_str, sep="    "))
    }
    extra <- (300000 - current_range) / 2
    new_start <- start_bp - extra
    new_end <- end_bp + extra
    wide_start <- new_start - 300000
    wide_end <- new_end + 300000
    thin_str <- paste("thin:", new_start, new_end, sep=" ")
    wide_str <- paste("wide:", wide_start, wide_end, sep=" ")
    return(paste(thin_str, wide_str, sep="    "))
}

get_bp_range(6940656, 6955947)



####### inspect indiv chromosomes  ----

rud_h_sigma150_pi_plot
# rud_both

rud_h_chr <- subset(rud_h_sigma150, Chr == "HG995199.1")
outliers0.00001_chr <- subset(outliers_0.00001, Chr == "HG995199.1")
outliers0.00001_chr <- outliers0.00001_chr %>%
    mutate(ymin=0, ymax=0.009,
           xmin = BP - sigma,
           xmax = BP + sigma) %>%
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>%
    select(-ymin_ymax)

outliers_0.005_chr <- outliers_0.005 %>%
    mutate(ymin=0, ymax=0.009,
           xmin = BP - sigma,
           xmax = BP + sigma) %>%
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>%
    select(-ymin_ymax)
low_outlier_0.005_chr <- subset(outliers_0.005_chr, Smoothed.Pi < bottom)


rud_h_chr_plot <- ggplot(rud_h_chr, aes(x=BP, y=Pi)) +
    
    # geom_vline(xintercept = 637182, col="red") +
    # geom_vline(xintercept = 806727, col="red") +
    # geom_vline(xintercept = 721954, col="red", alpha=0.5, size=3) +
    
    # outlier windows, pval < 0.00001
    geom_ribbon(data=outliers0.00001_chr,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="cyan",
                alpha=0.3) +

    geom_line(aes(group=Chr, x=BP, y=Smoothed.Pi), col="grey15") +
    geom_point(aes(x=BP, y=Smoothed.Pi), size=0.5) +
    
    # Custom the theme:
    theme_classic() +

    scale_y_continuous(expand = c(0, 0)) +
    
    labs(y=expression(pi),
         x="Position") +

    geom_vline(xintercept = 10686089, col="red") +
    geom_vline(xintercept = 5735417, col="red") +
    # geom_vline(xintercept = 1459310, col="red") +
    
    theme(legend.position="none",
          axis.text.x = element_text(),
          # axis.title.x = element_blank(),
          axis.text=element_text(size=10),
          axis.title=element_text(size=16))

rud_h_chr_plot


# based on genomic regions found above - blast to terr ref genome - get flat file - contains genes within region

## get genes from biobank flat files ----
# do separately for each file - so that know which genes from which region
rud_genes_flat <- data.frame(system("bash get_genes_from_flat.sh", intern=TRUE))
names(rud_genes_flat) <- "gene"
rud_genes_flat$gene <- lapply(strsplit(rud_genes_flat$gene, '"'), "[[", 2)
rud_genes_flat


## full genes list + gene names ----
rud_h_genes <- read.csv("../data/populations_smoothed/flat_genes_rud_h/smooth_pi_rud_h_genes.csv",
                        header=T)
rud_h_genes$id  <- 1:nrow(rud_h_genes)  # to keep orig order after merge

# gene name info
rud_h_gene_names <- read.csv("../data/populations_smoothed/flat_genes_rud_h/rud_h_gene_names_DAVID.tsv",
                              header=T, sep="\t")
rud_h_gene_names$Species <- NULL
names(rud_h_gene_names) <- c("gene", "name")

rud_h_genes <- merge(rud_h_genes, rud_h_gene_names, by="gene")
rud_h_genes <- rud_h_genes[order(rud_h_genes$id),]

# reorder columns
rud_h_genes <- rud_h_genes[,c("id", "chr_hort","bp_low_h", "bp_high_h", "chr_terr", "bp_low_t", "bp_high_t", "gene", "name")]

write.csv(rud_h_genes,
          file="../data/populations_smoothed/flat_genes_rud_h/smooth_pi_rud_h_genes_with_names.csv",
          row.names=F)


# length(intersect(GENES_IN_TERR$gene, rud_h_genes$gene))
# length(intersect(hort_h_gene_names$gene, rud_h_genes$gene))



## Fst ----

rud_h_fst_Os_Hd <- read.delim("../data/populations_smoothed/rud_h_smooth_sigma150/populations.fst_Os-Hd.tsv")
rud_h_fst_Ut_Hd <- read.delim("../data/populations_smoothed/rud_h_smooth_sigma150/populations.fst_Ut-Hd.tsv")
rud_h_fst_Ut_Os <- read.delim("../data/populations_smoothed/rud_h_smooth_sigma150/populations.fst_Ut-Os.tsv")

rud_h_fst_Os_Hd$compare <- c("Os-Hd")
rud_h_fst_Ut_Hd$compare <- c("Ut-Hd")
rud_h_fst_Ut_Os$compare <- c("Ut-Os")

rud_h_fst <- rbind(rud_h_fst_Os_Hd, rud_h_fst_Ut_Hd, rud_h_fst_Ut_Os)
rud_h_fst$compare <- as.factor(rud_h_fst$compare)

unique(rud_h_fst$Chr)
rud_h_fst$Chr <- as.factor(rud_h_fst$Chr)
rud_h_fst$Chr <- factor(rud_h_fst$Chr, levels=c(hort_linkage_groups))

## remove short unplaced linkage groups
rud_h_fst <- subset(rud_h_fst, Chr != "CAJOSO010000002.1" &
                        Chr != "CAJOSO010000007.1" &
                        Chr != "CAJOSO010000008.1")

# get cumulative basepair positions for plotting whole genome
rud_h_fst <- rud_h_fst %>% 
    
    # Compute chromosome size
    group_by(Chr) %>% 
    summarise(chr_len = max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(rud_h_fst, ., by=c("Chr"="Chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chr, BP) %>%
    mutate(BPcumul = BP + tot)

# get names of each chromosome (linkage group) and center of their bp for plotting x axis
axisdf = rud_h_fst %>% group_by(Chr) %>% summarize(center=( max(BPcumul) + min(BPcumul) ) / 2 )

# create dataframe for plotting striped background to separate linkage groups
strip_data_fst <- rud_h_fst %>% 
    select(BPcumul, Chr) %>% 
    group_by(Chr) %>% 
    mutate(ymin = 0, ymax = 0.12,
           xmin = min(BPcumul), xmax=max(BPcumul)) %>% 
    select(-BPcumul) %>% 
    distinct()
strip_data_fst$fill <- rep(c("a", "b"), length.out=nrow(strip_data_fst))

strip_data_fst <- strip_data_fst %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax)

# outliers
outliers_0.005 <- subset(rud_h_fst, Smoothed.AMOVA.Fst.P.value < 0.005)
outliers_0.001 <- subset(rud_h_fst, Smoothed.AMOVA.Fst.P.value < 0.001)
outliers_0.0005 <- subset(rud_h_fst, Smoothed.AMOVA.Fst.P.value < 0.0005)
outliers_0.0001 <- subset(rud_h_fst, Smoothed.AMOVA.Fst.P.value < 0.0001)
outliers_0.00001 <- subset(rud_h_fst, Smoothed.AMOVA.Fst.P.value < 0.00001)

length(unique(outliers_0.0001$X..Locus.ID))

outliers0.00001_unique <- unique(outliers_0.00001[c("X..Locus.ID", "compare")],)
table(outliers0.00001_unique$compare)


sigma=150000
outliers_0.00001_pm3sigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.12,
           xmin = BPcumul - sigma*3,
           xmax = BPcumul + sigma*3) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax)

outliers_0.00001_pm2sigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.12,
           xmin = BPcumul - sigma*2,
           xmax = BPcumul + sigma*2) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax)

outliers_0.00001_pmsigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.12,
           xmin = BPcumul - sigma,
           xmax = BPcumul + sigma) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax)


# plot - fst for each pop comparison in diff colour
rud_h_fst_plot <- ggplot(data=rud_h_fst, aes(x=BPcumul, y=Smoothed.AMOVA.Fst, col=compare)) +
    
    geom_ribbon(data=strip_data_fst,
                aes(y=y, xmin=xmin, xmax=xmax, group=Chr, fill=fill),
                inherit.aes=FALSE) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$Chr, breaks= axisdf$center,
                        expand=c(0.01,0.01)) +
    scale_y_continuous(expand = c(0, 0)) + 
    scale_fill_manual(name=NULL,
                      breaks=c("a", "b"),
                      # values=c("white", "#DFEAF9")) +     fcffa4                 values=c("white", "#DFEAF9")) +
                      values=c("white", "grey93")) +  # #fcffa4
    
    guides(fill="none") +
    
    geom_line(size=0.4) +

    scale_color_manual(name = expression(paste("Pairwise ", F["ST"])),
                       values = c("Os-Hd" = "#440154", 
                                  "Ut-Hd" = "#1f9e89",
                                  "Ut-Os" = "#b5de2b"),
                       labels = c("Os-Hd" = "Ouse - Hillesden", 
                                  "Ut-Hd" = "Upton - Hillesden",
                                  "Ut-Os" = "Upton - Ouse")) + 

    annotate(geom="text", x=c(outliers_0.00001$BPcumul), y=0.115, label="↓", 
             col="red", size=8) +
    
    theme_classic() +
    labs(y=expression(F["ST"]),
         x=expression(paste("\nBasepair position along ", italic("B. hortorum"), " linkage groups"))) +
    theme(
        # legend.position="none",
          legend.title = element_text(size=25),
          legend.text = element_text(size=25),
          axis.text.x = element_text(angle = 90), 
          axis.text=element_text(size=20),
          axis.title=element_text(size=25),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          legend.position = "bottom",
          plot.margin = margin(t = 12, r = 5.5, b = 5.5, l = 5.5, unit = "pt")) +
    guides(color = guide_legend(override.aes = list(size = 8)))

rud_h_fst_plot


grDevices::cairo_pdf("../data/populations_smoothed/rud_h_fst.pdf",
                     width=16)
rud_h_fst_plot
dev.off()

# ggsave(plot=rud_h_fst_plot,
#        file="../data/populations_smoothed/rud_h_fst.pdf",
#        width=16)

# rud_h_sigma150_plots <- grid.arrange(rud_h_pi_plot_for_combo, rud_h_fst_plot, heights=c(1,1.3))

# library(cowplot)
# library(grDevices)

rud_both <- ggdraw() +
    draw_plot(rud_h_pi_plot_for_combo, x = 0, y = 0.6, width = 1, height = 0.4) +
    draw_plot(rud_h_fst_plot, x = 0, y = 0, width = 1, height = 0.6) +
    draw_plot_label(label = c("A", "B"), size = 25, x = c(0,0), y = c(1,0.6))
rud_both

grDevices::cairo_pdf("../thesis/figures/rud_pi_fst.pdf",
                     height=15, width=16)
rud_both
dev.off()

grDevices::cairo_pdf("../data/populations_smoothed/rud_h_pi_fst.pdf",
                     height=15, width=16)
rud_both
dev.off()



## inspect indiv fst chr ----

rud_h_fst_plot
# rud_both

rud_h_chr <- subset(rud_h_fst, Chr == "HG995188.1")
outliers0.00001_chr <- subset(outliers_0.00001, Chr == "HG995188.1")
outliers0.00001_chr <- outliers0.00001_chr %>%
    mutate(ymin=0, ymax=0.1,
           xmin = BP - sigma,
           xmax = BP + sigma) %>%
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>%
    select(-ymin_ymax)

rud_h_chr_plot <- ggplot(data=rud_h_chr, aes(x=BP, y=Smoothed.AMOVA.Fst, col=compare)) +
    
    scale_y_continuous(expand = c(0, 0)) + 
    geom_line(size=0.4) +

    scale_color_manual(name = expression(paste("Pairwise ", F["ST"])),
                   values = c("Os-Hd" = "#440154", 
                              "Ut-Hd" = "#1f9e89",
                              "Ut-Os" = "#b5de2b"),
                   labels = c("Os-Hd" = "Ouse - Hillesden", 
                              "Ut-Hd" = "Upton - Hillesden",
                              "Ut-Os" = "Upton - Ouse")) + 
    
    annotate(geom="text", x=c(outliers0.00001_chr$BP), y=0.115, label="↓", 
             col="red", size=8) +
    
    # ## bayescan
    geom_vline(xintercept = 20927101, col="red") +
    # geom_vline(xintercept = 3251256, col="red") +
    # geom_vline(xintercept = 20927101, col="red") +
    # 
    # ## pcadapt
    # geom_vline(xintercept = 402530, col="black") +
    # geom_vline(xintercept = 426308, col="black") +
    # geom_vline(xintercept = 9494806, col="black") +
    # geom_vline(xintercept = 14951165, col="black") +
    
    theme_classic()

rud_h_chr_plot

