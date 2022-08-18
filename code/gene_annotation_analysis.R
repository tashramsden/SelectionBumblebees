## analysis of gene ontology ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(scales)
library(VennDiagram)
library(grDevices)
library(ggsignif)


## pi compare between spp ----
rud_pi <- rud_h_sigma150
rud_pi$spp <- "rud"
hort_pi <- hort_h_sigma150
hort_pi$spp <- "hort"
terr_pi <- terr_t_sigma150
terr_pi$spp <- "terr"

all_pi <- rbind(rud_pi, hort_pi, terr_pi)

mod <- aov(all_pi$Smoothed.Pi ~ all_pi$spp)
summary(mod)
TukeyHSD(mod)
plot(TukeyHSD(mod, conf.level=.95), las = 2)

# ANOVA: spp is sig factor in nucleotide diversity
# Tukey: all 3 spp pi sig diff - B. ruderatus lower than other 2

all_pi$spp <- as.factor(all_pi$spp)
all_pi$spp <- factor(all_pi$spp, levels=c("rud", "hort", "terr"))

# make separate dataframe of outliers to add jitter to plot
outliers <- 
    all_pi %>%
    group_by(spp) %>%
    filter(Smoothed.Pi > quantile(Smoothed.Pi, 0.75) + 1.5 * IQR(Smoothed.Pi) | 
               Smoothed.Pi < quantile(Smoothed.Pi, 0.25) - 1.5 * IQR(Smoothed.Pi))

nuc_div_spp_plot <- ggplot(all_pi, aes(x=spp, y=Smoothed.Pi, col=spp)) +
    geom_boxplot(alpha=0.05, outlier.shape=NA, size=1.5) +
    geom_jitter(height = 0, width = 0.2, data = outliers, alpha=0.3) +
    scale_y_continuous(expand=c(0,0.0008)) +
    scale_x_discrete(labels=c(expression(italic("B. ruderatus")),
             expression(italic("B. hortorum")),
             expression(italic("B. terrestris")))) +
    labs(x="Species", y="Nucleotide Diversity") +
    theme_classic() +
    theme(axis.title = element_text(size=40),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin=margin(r=15)),
          axis.text = element_text(size=35),
          legend.position = "none") +
    geom_signif(comparisons = list(c("rud", "hort"),
                                   c("rud", "terr"),
                                   c("hort", "terr")),
                map_signif_level=TRUE, textsize = 12, col="black")
nuc_div_spp_plot

ggsave(nuc_div_spp_plot,
       file="../presentation/nuc_div_spp.jpg",
       width=25, height=15)

nuc_div_spp_plot <- nuc_div_spp_plot + 
    theme(axis.title = element_text(size=25),
          axis.title.x = element_text(margin=margin(t=10)),
          axis.title.y = element_text(margin=margin(r=10)),
          axis.text = element_text(size=20))
nuc_div_spp_plot
ggsave(nuc_div_spp_plot,
       file="../thesis/figures/supplementary_nuc_div_spp.pdf",
       width=15, height=10)


## genes found in regions of signif elev smoothed pi ----

terr_genes <- read.csv("../data/populations_smoothed/flat_genes_terr_t/smooth_pi_terr_t_genes_with_names_both.csv",
                       header=T)

rud_genes <- read.csv("../data/populations_smoothed/flat_genes_rud_h/smooth_pi_rud_h_genes_with_names.csv",
                      header=T)

hort_genes <- read.csv("../data/populations_smoothed/flat_genes_hort_h/smooth_pi_hort_h_genes.csv",
                       header=T)

all_genes <- c(terr_genes$gene, hort_genes$gene, rud_genes$gene)
length(unique(all_genes))


# venn diagram of genes in each spp
myCol <- c("#fde725", "#5ec962", "#3b528b")  # viridis colours
temp <- venn.diagram(
    x = list(terr_genes$gene, hort_genes$gene, rud_genes$gene),
    category.names = c("B. terrestris", "B. hortorum", "B. ruderatus"),
    # circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    # numbers
    cex = 3,
    fontfamily = "sans",
    # categories
    cat.cex = 3,
    cat.fontfamily = "sans",
    cat.fontface = "italic",
    
    cat.default.pos = "outer",
    cat.pos = c(-26, 26, 180),

    filename = NULL
)

pdf(file="../thesis/figures/venn_pi.pdf", width=8, height=8)
grid.draw(temp)
dev.off()


# intersection between all 3 species:
length(Reduce(intersect, list(rud_genes$gene, hort_genes$gene, terr_genes$gene)))
# 15 genes found in all 3 species:
Reduce(intersect, list(rud_genes$gene, hort_genes$gene, terr_genes$gene))

# rud + hort intersect:
rud_hort <- intersect(rud_genes$gene, hort_genes$gene)  # 21 shared total
setdiff(rud_hort, terr_genes$gene)  # 6 genes shared hort + rud and not terr

# rud + terr intersect:
rud_terr <- intersect(rud_genes$gene, terr_genes$gene)  # 17 shared genes
setdiff(rud_terr, hort_genes$gene)  # 2 genes shared rud + terr not hort

# hort + terr intersect:
hort_terr <- intersect(hort_genes$gene, terr_genes$gene)  # 20 shared genes
setdiff(hort_terr, rud_genes$gene)  # 5 genes shared hort + terr not rud


# Functional categories - pie charts ----
pie_cols <- c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fca50a")

## B. ruderatus functional categories ----

rud_cat <- read.csv("../data/populations_smoothed/flat_genes_rud_h/rud_h_genes_categories.csv",
                    header=T)

ggplot(rud_cat, aes(x="", fill=category)) +
    geom_bar(width=1, color="white") +
    theme_classic() +
    scale_fill_viridis(discrete = TRUE) +
    coord_polar("y", start=0)

table(rud_cat$category)

rud_cat <- rud_cat %>% 
    mutate(cat_simple = case_when(category == "uncharacterised" ~ "Uncharacterised",
                                  category == "lncRNA" ~ "Long non-coding RNA",
                                  category == "immunity" ~ "Immunity",
                                  category == "morphogenesis" ~ "Morphogenesis",
                                  category == "neurology" ~ "Neurology",
                                  TRUE ~ "Other"))

rud_cat_simple <- data.frame(table(rud_cat$cat_simple))
names(rud_cat_simple) <- c("func_cat", "count")


rud_cat_pie <- ggplot(rud_cat_simple, aes(x="", y=count, fill=func_cat)) +
    geom_bar(stat="identity", width=1, color="white") +
    theme_classic() +
    scale_fill_manual(values=pie_cols, name="Functional Category") + 
    # scale_fill_viridis(discrete = TRUE, name="Functional Category", option = "viridis") +
    labs(x="", y="Count", title=expression(italic("B. ruderatus"))) +
    geom_text(aes(label = count, x=1.3),
              position = position_stack(vjust = 0.5),
              size=10, col="white") +
    theme_void() + # remove background, grid, numeric labels 
    coord_polar("y", start=0) +
    theme(title = element_text(size=25),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t = 12, r = 0, b = 0, l = 5.5, unit = "pt"))
          # legend.position = "none")
rud_cat_pie



## B. terr function categories ----

terr_cat <- read.csv("../data/populations_smoothed/flat_genes_terr_t/terr_t_genes_categories.csv",
                    header=T, skip=3)

ggplot(terr_cat, aes(x="", fill=category)) +
    geom_bar(width=1, color="white") +
    theme_classic() +
    scale_fill_viridis(discrete = TRUE) +
    coord_polar("y", start=0)

table(terr_cat$category)

terr_cat <- terr_cat %>% 
    mutate(cat_simple = case_when(category == "uncharacterised" ~ "Uncharacterised",
                                  category == "lncRNA" ~ "Long non-coding RNA",
                                  category == "immunity" ~ "Immunity",
                                  category == "morphogenesis" ~ "Morphogenesis",
                                  category == "neurology" ~ "Neurology",
                                  TRUE ~ "Other"))

terr_cat_simple <- data.frame(table(terr_cat$cat_simple))
names(terr_cat_simple) <- c("func_cat", "count")

terr_cat_pie <- ggplot(terr_cat_simple, aes(x="", y=count, fill=func_cat)) +
    geom_bar(stat="identity", width=1, color="white") +
    theme_classic() +
    scale_fill_manual(values=pie_cols, name="Functional Category") +
    # scale_fill_viridis(discrete = TRUE, name="Functional Category", option = "viridis") +
    labs(x="", y="Count", title=expression(italic("B. terrestris"))) +
    geom_text(aes(label = count, x=1.3),
              position = position_stack(vjust = 0.5),
              size=10, col="white") +
    theme_void() + # remove background, grid, numeric labels 
    coord_polar("y", start=0) +
    theme(title = element_text(size=25),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t = 12, r = 5.5, b = 0, l = 5.5, unit = "pt"))

# legend.position = "none")
terr_cat_pie



## B. hort function categories ----

hort_cat <- read.csv("../data/populations_smoothed/flat_genes_hort_h/hort_h_genes_categories.csv",
                     header=T, skip=3)

ggplot(hort_cat, aes(x="", fill=category)) +
    geom_bar(width=1, color="white") +
    theme_classic() +
    scale_fill_viridis(discrete = TRUE) +
    coord_polar("y", start=0)

table(hort_cat$category)

hort_cat <- hort_cat %>% 
    mutate(cat_simple = case_when(category == "uncharacterised" ~ "Uncharacterised",
                                  category == "lncRNA" ~ "Long non-coding RNA",
                                  category == "immunity" ~ "Immunity",
                                  category == "morphogenesis" ~ "Morphogenesis",
                                  category == "neurology" ~ "Neurology",
                                  TRUE ~ "Other"))

hort_cat_simple <- data.frame(table(hort_cat$cat_simple))
names(hort_cat_simple) <- c("func_cat", "count")

hort_cat_pie <- ggplot(hort_cat_simple, aes(x="", y=count, fill=func_cat)) +
    geom_bar(stat="identity", width=1, color="white") +
    theme_classic() +
    scale_fill_manual(values=pie_cols, name="Functional Category") +
    # scale_fill_viridis(discrete = TRUE, name="Functional Category", option = "viridis") +
    # scale_fill_brewer(palette="Set1", name="Functional Category") +
    labs(x="", y="Count", title=expression(italic("B. hortorum"))) +
    geom_text(aes(label = count, x=1.3),
              position = position_stack(vjust = 0.5),
              size=10, col="white") +
    theme_void() + # remove background, grid, numeric labels 
    coord_polar("y", start=0) +
    theme(title = element_text(size=25),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t = 12, r = 5.5, b = 0, l = 5.5, unit = "pt"))
# legend.position = "none")
hort_cat_pie


# all pies ----
pies <- plot_grid( terr_cat_pie + theme(legend.position="none"),
                   hort_cat_pie + theme(legend.position="none"),
                   rud_cat_pie + theme(legend.position="none"),
                   align = 'vh',
                   nrow = 1
)
pies

# all have same legend
pie_legend <- get_legend(terr_cat_pie + theme(legend.position="bottom",
                                              legend.title = element_text(size=30),
                                              legend.text = element_text(size=30)))
                         
# all pies w legend at bottom
cat_pies2 <- ggdraw() +
    draw_plot(pies, x = 0, y = 0.1, width = 1, height = 0.9) +
    draw_plot(pie_legend, x = 0, y = 0, width = 1, height = 0.3)
    # draw_plot_label(label = c("A", "B"), size = 25, x = c(0,0), y = c(1,0.6))
cat_pies2

ggsave(cat_pies2,
       file = "../thesis/figures/cat_pies.pdf",
       width=16, height=10)

# ggsave(cat_pies2,
#        file = "../presentation/cat_pies.jpg",
#        width=16, height=10)


## genes found in colgan et al 2019 and bebane 2019 to be diff expressed on neonicotinoid exposure ----

bebane_pests_actual <- read_excel("../other/bebane_et_al_2019.xlsx",  # is same as above!
                           sheet = "differentially_expressed_genes")
intersect(bebane_pests_actual$LOC, terr_genes$gene) # "LOC100643070"-ryanodine receptor-immunity 
                                                    # "LOC100631088"-phosphoenolpyruvate carboxykinase 
                                                    # "LOC100650142"-furin-like protease 1-immuntiy
intersect(bebane_pests_actual$LOC, hort_genes$gene) # "LOC100644174"-flocculation protein FLO11-immunity 
                                                    # "LOC100646870"-unchar 
                                                    # "LOC100644143"-putative inorganic phosphate cotransporter 
                                                    # "LOC100646748"-GTP-binding nuclear protein Ran 
                                                    # "LOC100648817"-protein BUD31 homolog
# "LOC
intersect(bebane_pests_actual$LOC, rud_genes$gene)  # "LOC100650874"-unchar 
                                                    # "LOC100644830"-cleavage stimulation factor subunit 1 - pre-mRNA formation

# colgan et al 2019 (not 2022) diff expressed
colgan_2019_pests <- read_excel("../other/colgan_et_al_supplementary/Table_S9_pesticide_response_genes.xlsx",
                           sheet = "Table_S9c_colgan_worker_DE",
                           skip = 7)
colgan_2019_pests <- colgan_2019_pests[,1:2]

intersect(colgan_2019_pests$`NCBI Gene symbol`, terr_genes$gene) # "LOC100631088"
intersect(colgan_2019_pests$`NCBI Gene symbol`, hort_genes$gene) # none
intersect(colgan_2019_pests$`NCBI Gene symbol`, rud_genes$gene) # "LOC100651821"-latrophilin-like protein LAT-2 - neurology
