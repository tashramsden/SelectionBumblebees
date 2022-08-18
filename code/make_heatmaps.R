# heat maps! ----

rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(ggpubr)
library(cowplot)

## ruderatus ----

rud_haps_outliers <- read.table("../data/bayescan/rud_sites/outliers_alleles_sites.csv",
                                header=T)

# find dominant allele freq at each locus
dom_freqs <- rud_haps_outliers %>% 
    group_by(site, locus_id) %>% 
    summarise(most_freq = tail(names(sort(table(alleles))), 1),
              dom_freq = sum(alleles == most_freq) / length(alleles))

unique(dom_freqs$locus_id)

# loci in order of qval (lowest to highest - ie most significant down)
dom_freqs$locus_id <- as.factor(dom_freqs$locus_id)
dom_freqs$locus_id <- factor(dom_freqs$locus_id,
                             levels=c("64429","1702","64341","1709","22012","62734","5335"))

# geo locations west to east
dom_freqs$site <- as.factor(dom_freqs$site)
dom_freqs$site <- factor(dom_freqs$site,
                         levels=c("Upton", "Hillesden", "Ouse"))

# add gene annotation info
bayescan_results <- read.csv("../data/bayescan/rud_bayescan_genes_cats.csv")
bayescan_results <- bayescan_results[, 1:9]
names(bayescan_results)[names(bayescan_results) == "hort_aligned_locus"] <- "locus_id"
locus_gene <- bayescan_results[,c("locus_id", "gene")]
locus_gene$locus_id <- as.factor(locus_gene$locus_id)

dom_freqs <- dom_freqs %>% 
    left_join(locus_gene, by="locus_id")

bayes_loci <- unique(dom_freqs$locus_id)

# loci in order of qval (lowest to highest - ie most significant down)
dom_freqs$gene <- as.factor(dom_freqs$gene)
dom_freqs$gene <- factor(dom_freqs$gene,
                             levels=c("LOC100651177",
                                      "No gene 1",
                                      "LOC100647201",
                                      "No gene 2",
                                      "LOC100643645",
                                      "LOC100649273",
                                      "LOC100643486"))

# only want snps detected by >1 pop-diff test
dom_freqs <- subset(dom_freqs, gene != "LOC100643486")

# plot heatmap of dom allele freq at each locus and geo location
rud_heatmap_genes <- ggplot(dom_freqs, 
                            aes(x=site, y=reorder(gene, desc(gene)),
                                fill=dom_freq)) +
    geom_tile(color="white", size=0.3) +
    scale_x_discrete(expand=c(0,0), position = "top") +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_viridis(discrete=FALSE, direction=-1, option="viridis",
                       name="Major allele \n frequency", limits=c(0.5, 1)) +
    # scale_fill_gradient(low = "white", high = "blue") +
    coord_fixed(ratio=0.4) +
    labs(x="", y="Genes containing outlier SNPs", ) +
    theme_classic() +
    # annotation_custom(tableGrob(tests_per_gene1, rows=NULL, theme=mytheme),
    # xmin=6, xmax=6, ymin=1, ymax=8) +
    theme(axis.title = element_text(size=25),
          axis.title.y = element_text(margin=margin(t = 0, r = 15, b = 0, l = 0)),
          axis.text = element_text(size=25),
          # axis.text.x = element_text(vjust=-10),
          legend.title = element_text(size=15),
          legend.text = element_text(size=15),)
rud_heatmap_genes

# remove legend for combination plot
rud_heatmap_no_legend <- rud_heatmap_genes + theme(legend.position = "none")

# dataframe for plotting of which loci were found by which pop-diff selection tests
tests_per_gene <- data.frame(
                            # Gene = c("LOC100651177", 
                            #           "No gene 1",    
                            #           "LOC100647201", 
                            #           "No gene 2",    
                            #           "LOC100643645", 
                            #           "LOC100649273", 
                            #           "LOC100643486"),
                             SmoothFst = c("X",
                                              "X",
                                              "X",
                                              "X",
                                              "X",
                                              ""),
                                              # ""),
                             BayeScan = c("X",
                                          "X",
                                          "X",
                                          "X",
                                          "X",
                                          "X"),
                                          # "X"),
                             PCAdapt = c("X",
                                         "",
                                         "X",
                                         "",
                                         "",
                                         "X"))
                                         # ""))

names(tests_per_gene) <- c("Smooth Fst", "BayeScan", "PCAdapt")

# make pretty table of tests per locus
heat_tests <- ggtexttable(tests_per_gene, rows = NULL, 
                        theme = ttheme("blank",  
                                       colnames.style = colnames_style(color = "grey40", fill="white",
                                                                       size=24, face="plain",
                                                                       vjust=1.2),
                                       base_size = 25,
                                       padding = unit(c(10, 12.5), "mm"))
                        ) %>% 
    table_cell_font(row=2:7, column=1:3, color = "grey30", size=20) %>%
    tab_add_hline(at.row = 2:7, row.side = "top", linewidth = 2, linetype=2) %>% 
    tab_add_vline(at.column = 1:2, column.side = "right", from.row=2, linewidth = 3, linetype=2)

heat_tests

# combine heatmap and table of tests
heats2 <- ggdraw() +
    draw_plot(rud_heatmap_no_legend, x = -0.05, y = 0.1, width = 0.7, height = 0.9) +
    draw_plot(heat_tests, x = 0.66, y = 0.048, width = 0.2, height = 1)
heats2

# ggsave(heats2,
#        file="../thesis/figures/rud_heatmap.pdf",
#        width=16)

# add legend back to combination plot
heat_legend <- get_legend(rud_heatmap_genes + theme(legend.position="right",
                                                    legend.title = element_text(size=25),
                                                    legend.text = element_text(size=25),
                                                    legend.key.size = unit(1.2, 'cm')))

heat_full <- ggdraw() +
    draw_plot(heats2, x = 0, y = 0.0, width = 0.9, height = 1) +
    draw_plot(heat_legend, x = 0.78, y = -0.45, width = 0.3, height = 2)
heat_full

ggsave(heat_full,
       file="../thesis/figures/rud_heatmap.pdf",
       width=18)

# compare dominant allele freq at the geo locations
dom_freqs %>% group_by(site) %>% 
    summarise(mean=mean(dom_freq),
              se=sd(dom_freq) / sqrt(length(dom_freq)))

# anova + tukey to test significance
mod <- aov(dom_freqs$dom_freq ~ dom_freqs$site)
summary(mod)
TukeyHSD(mod)
plot(TukeyHSD(mod, conf.level=.95), las = 2)
