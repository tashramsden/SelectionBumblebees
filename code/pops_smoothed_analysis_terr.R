# stacks smoothed stats + bootstrapping evaluation - B. terrestris

rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(readxl)
library(scales)
library(cowplot)


## NUCLEOTIDE DIVERSITY - PI ----

terr_linkage_groups_t <- c("NC_063269.1",
                         "NC_063270.1",
                         "NC_063271.1",
                         "NC_063272.1",
                         "NC_063273.1",
                         "NC_063274.1",
                         "NC_063275.1",
                         "NC_063276.1",
                         "NC_063277.1",
                         "NC_063278.1",
                         "NC_063279.1",
                         "NC_063280.1",
                         "NC_063281.1",
                         "NC_063282.1",
                         "NC_063283.1",
                         "NC_063284.1",
                         "NC_063285.1",
                         "NC_063286.1",
                         "NW_025963548.1",  # unplaced
                         "NW_025963549.1",
                         "NW_025963552.1",
                         "NW_025963557.1",
                         "NW_025963558.1",
                         "NW_025963586.1",
                         "NW_025963599.1",
                         "NW_025963600.1",
                         "NW_025963605.1")


## B. terrestris - aligned to terr ref genome ----

## sigma150 ----
terr_t_sigma150 <- read.delim("../data/populations_smoothed/terr_t_smooth_sigma150/populations.sumstats.tsv",
                             skip=4)

unique(terr_t_sigma150$Chr)
terr_t_sigma150$Chr <- as.factor(terr_t_sigma150$Chr)
terr_t_sigma150$Chr <- factor(terr_t_sigma150$Chr, levels=c(terr_linkage_groups_t))
# levels(terr_t_sigma150$Chr)

# remove data where smoothed pi = -1 (here overlapping loci not incl in calculation and given val -1)
terr_t_sigma150 <- subset(terr_t_sigma150, Smoothed.Pi != -1)

## remove short unplaced linkage groups
terr_t_sigma150 <- subset(terr_t_sigma150, Chr != "NW_025963548.1" &
                              Chr != "NW_025963549.1" &
                              Chr != "NW_025963552.1" &
                              Chr != "NW_025963557.1" &
                              Chr != "NW_025963558.1" &
                              Chr != "NW_025963586.1" &
                              Chr != "NW_025963599.1" &
                              Chr != "NW_025963600.1" &
                              Chr != "NW_025963605.1")

# get cumulative basepair positions for plotting whole genome
terr_t_sigma150 <- terr_t_sigma150 %>% 
    
    # Compute chromosome size
    group_by(Chr) %>% 
    summarise(chr_len = max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(terr_t_sigma150, ., by=c("Chr"="Chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chr, BP) %>%
    mutate(BPcumul = BP + tot)

# get names of each chromosome (linkage group) and center of their bp for plotting x axis
axisdf = terr_t_sigma150 %>% group_by(Chr) %>% summarize(center=( max(BPcumul) + min(BPcumul) ) / 2 )

outliers_0.05 <- subset(terr_t_sigma150, Smoothed.Pi.P.value < 0.05)
outliers_0.01 <- subset(terr_t_sigma150, Smoothed.Pi.P.value < 0.01)
outliers_0.005 <- subset(terr_t_sigma150, Smoothed.Pi.P.value < 0.005)
outliers_0.001 <- subset(terr_t_sigma150, Smoothed.Pi.P.value < 0.001)
outliers_0.0005 <- subset(terr_t_sigma150, Smoothed.Pi.P.value < 0.0005)
outliers_0.0001 <- subset(terr_t_sigma150, Smoothed.Pi.P.value < 0.0001)
outliers_0.00001 <- subset(terr_t_sigma150, Smoothed.Pi.P.value < 0.00001)

length(unique(outliers_0.00001$X..Locus.ID))

write.table(outliers_0.005, file="../data/populations_smoothed/terr_t_smooth_sigma150/outliers_0.005")
write.table(outliers_0.00001, file="../data/populations_smoothed/terr_t_smooth_sigma150/outliers_0.00001")

# hist(terr_t_sigma150$Smoothed.Pi.P.value)
hist(terr_t_sigma150$Smoothed.Pi)
min(outliers_0.05$Smoothed.Pi)

top <- quantile(terr_t_sigma150$Smoothed.Pi, 0.95)
bottom <- quantile(terr_t_sigma150$Smoothed.Pi, 0.05)
low <- terr_t_sigma150[which(terr_t_sigma150$Smoothed.Pi < bottom),]

# create dataframe for plotting striped background to separate linkage groups
strip_data <- terr_t_sigma150 %>% 
    select(BPcumul, Chr) %>% 
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
    select(-ymin_ymax)  %>% 
    mutate(shared = case_when(Chr == "NC_063278.1" & BP <= 16607539 ~ "red",  # shared all spp
                              TRUE ~ "#fca50a"))  # not shared by all spp


outliers_0.00001_pm2sigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.009,
           xmin = BPcumul - sigma*2,
           xmax = BPcumul + sigma*2) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax)  %>% 
    mutate(shared = case_when(Chr == "NC_063278.1" & BP <= 16607539 ~ "red",  # shared all spp
                              TRUE ~ "#fca50a"))  # not shared by all spp


outliers_0.00001_pmsigmawindow <- outliers_0.00001 %>% 
    mutate(ymin=0, ymax=0.009,
           xmin = BPcumul - sigma,
           xmax = BPcumul + sigma) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax)  %>% 
    mutate(shared = case_when(Chr == "NC_063278.1" & BP <= 16607539 ~ "red",  # shared all spp
                              TRUE ~ "#fca50a"))  # not shared by all spp


terr_t_sigma150_pi_plot <- ggplot(terr_t_sigma150, aes(x=BPcumul, y=Pi)) +
    
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
                fill="#fca50a", # outliers_0.00001_pm3sigmawindow$shared,
                alpha=0.2) +
    geom_ribbon(data=outliers_0.00001_pm2sigmawindow,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="#fca50a") + # outliers_0.00001_pm2sigmawindow$shared) +
    geom_ribbon(data=outliers_0.00001_pmsigmawindow,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="#fca50a") + # outliers_0.00001_pmsigmawindow$shared) +
    
    geom_line(aes(group=Chr, BPcumul, Smoothed.Pi), col="grey15") +
    # geom_point(aes(group=Chr, BPcumul, Smoothed.Pi), col="grey15", size=0.25) +
    
    annotate(geom="point", x=183550000, y=0.0087, shape=25, col="blue", fill="blue", size=5) +
    
    # Custom the theme:
    theme_classic() +

    labs(y=expression(pi),
         x=expression(paste("Basepair position along ", italic("B. terrestris"), " linkage groups"))) +
         # expression(paste("No. of ", italic("bacteria X"), " isolates with corresponding types"))) +
    theme(legend.position="none",
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.text.x = element_text(angle = 90),
          # axis.title.x = element_blank(),
          axis.text=element_text(size=20),
          axis.title=element_text(size=25))
terr_t_sigma150_pi_plot


mean(terr_t_sigma150$Smoothed.Pi)
se = sd(terr_t_sigma150$Smoothed.Pi) / sqrt(length(terr_t_sigma150$Smoothed.Pi))
se


# terr_t_outliers0.005 <- read.table("../data/populations_smoothed/terr_t_smooth_sigma150/outliers_0.005")
# terr_t_outliers0.00001 <- read.table("../data/populations_smoothed/terr_t_smooth_sigma150/outliers_0.00001")

ggsave(plot=terr_t_sigma150_pi_plot,
       file="../data/populations_smoothed/terr_t_pi.pdf",
       width=16)

ggsave(plot=terr_t_sigma150_pi_plot,
       file="../thesis/figures/terr_pi.pdf",
       width=16, height=8)


## for presentation ----

terrplot <- terr_t_sigma150_pi_plot +
    # scale_x_discrete(label=seq(1, 18, 1))
    # theme(axis.text.x = element_text(c(seq(1, 18, 1))))
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
terrplot

ggsave(plot=terrplot,
       file="../presentation/terr_pi.jpg",
       width=18, height=4)




####### JUST CHR 1 TO INSPECT CONSERVED REGION COLGAN ET AL ----
terr_chr1 <- subset(terr_t_sigma150, Chr == "NC_063269.1") 
outliers0.00001_chr1 <- subset(outliers_0.00001, Chr == "NC_063269.1")
outliers0.00001_chr1 <- outliers0.00001_chr1 %>% 
    mutate(ymin=0, ymax=0.005,
           xmin = BPcumul - sigma,
           xmax = BPcumul + sigma) %>% 
    pivot_longer(cols=c(ymin, ymax), values_to = "y", names_to="ymin_ymax") %>% 
    select(-ymin_ymax)

terr_chr1 <- terr_chr1 %>% 
    mutate(BP_plot = BP / 1000000)


terr_chr1_plot <- ggplot(terr_chr1, aes(x=BP, y=Pi)) +

    geom_vline(xintercept = 637182, col="#7ad151") +
    geom_vline(xintercept = 806727, col="#7ad151") +
    geom_vline(xintercept = 721954, col="#7ad151", alpha=0.5, size=5) +
    
    # outlier windows, pval < 0.00001
    geom_ribbon(data=outliers0.00001_chr1,
                aes(y=y, xmin=xmin, xmax=xmax, group=X..Locus.ID),
                inherit.aes=FALSE,
                fill="#fca50a",
                alpha=0.3) +
    
    geom_line(aes(group=Chr, BPcumul, Smoothed.Pi), col="grey15") +
    geom_point(aes(x=BP, y=Smoothed.Pi), size=0.5) +
    
    # Custom the theme:
    theme_classic() +
    
    scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6),
                       breaks=seq(0, max(terr_chr1$BP), 1000000),
                       expand = c(0.01, 0.01)) + 
    scale_y_continuous(expand = c(0, 0)) +
    
    labs(y=expression(pi),
         x="Position (Mb)") +
    theme(legend.position="none",
          axis.text.x = element_text(),
          # axis.title.x = element_blank(),
          axis.text=element_text(size=20),
          axis.title=element_text(size=25),
          # add little bit padding at top to fit in axis text!
          plot.margin = margin(t = 12, r = 5.5, b = 5.5, l = 5.5, unit = "pt"))

terr_chr1_plot

ggsave(plot=terr_chr1_plot,
       file="../data/populations_smoothed/terr_chr1.pdf",
       width=16)

ggsave(plot=terr_chr1_plot,
       file="../thesis/figures/supplementary_terr_chr1.pdf",
       width=16)


terr_both <- ggdraw() +
    draw_plot(terr_t_sigma150_pi_plot, x = 0, y = 0.4, width = 1, height = 0.6) +
    draw_plot(terr_chr1_plot, x = 0, y = 0, width = 0.5, height = 0.4) +
    # draw_plot(p2, x = 0, y = 0, width = 1, height = 0.5) +
    draw_plot_label(label = c("A", "B: NC_063269.1"), size = 15, x = c(0,0), y = c(1,0.4))
terr_both

# ggsave(plot=terr_both,
#        file="../data/populations_smoothed/terr_arrange.pdf",
#        width=16)


## full genes list + gene names ----
terr_t_genes <- read.csv("../data/populations_smoothed/flat_genes_terr_t/smooth_pi_terr_t_genes_with_names.csv",
                        header=T)
terr_t_genes$id  <- 1:nrow(terr_t_genes)  # to keep orig order after merge

# gene name info
terr_t_gene_names <- read.csv("../data/populations_smoothed/flat_genes_terr_t/terr_t_gene_names_DAVID.tsv",
                             header=T, sep="\t")
terr_t_gene_names$Species <- NULL
names(terr_t_gene_names) <- c("gene", "DAVIDname")

terr_t_gene_names <- head(terr_t_gene_names, -3)

terr_t_genes <- merge(terr_t_genes, terr_t_gene_names, by="gene")
terr_t_genes <- terr_t_genes[order(terr_t_genes$id),]

# reorder columns
terr_t_genes <- terr_t_genes[,c("id", "chr", "bp_low", "bp_high", "gene", "name", "DAVIDname")]

write.csv(terr_t_genes,
          file="../data/populations_smoothed/flat_genes_terr_t/smooth_pi_terr_t_genes_with_names_both.csv",
          row.names=F)
