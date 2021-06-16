library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
library(tidyr)
library(ochRe)
library(ggpubr)
library(cowplot)
library(scales)

source(file = "../../config.txt")

###################
## Figure 3A
###################

# Count the number of genes in which each CAZY module is annotated

cazy <- read.delim(file = paste0(mydir, "/Module_4/cazy_05/hmmer.out"))
cazy %>% 
  select(Gene.ID, HMM.Profile) %>% 
  setNames(., c("gene", "HMM")) -> cazy
cazy$gene <- str_split(cazy$gene, "_", simplify = TRUE)[,1]
cazy %>% 
  group_by(gene) %>% 
  distinct(HMM, .keep_all = TRUE) -> cazy_unique_HMM
cazy_unique_HMM$HMM <- sub(".hmm", "", cazy_unique_HMM$HMM)
cazy_unique_HMM$family <- gsub('[[:digit:]]+', '', cazy_unique_HMM$HMM)
cazy_unique_HMM$family <- str_split(cazy_unique_HMM$family, "_", simplify = TRUE)[,1]

# Add info about up / or downregulation

gene_expression <- read.table(file = paste0(mydir, "/Tables/gene.expression.out"))
gene_expression %>% 
  mutate(up = ifelse(gene_expression$log2FoldChange_as > 1 & gene_expression$padj_as < 0.05,
                     TRUE,
                     FALSE)) %>% 
  mutate(down = ifelse(gene_expression$log2FoldChange_as < -1 & gene_expression$padj_as < 0.05,
                     TRUE,
                     FALSE)) %>% 
  select(gene, up, down) -> gene_expression_cats

# Count up and downregulated genes
sum(gene_expression_cats$up)
sum(gene_expression_cats$down)

cazy_unique_HMM %>% 
  left_join(gene_expression_cats) -> cazy_unique_HMM_expression

# Omit cazys for which there is no gene expression data 
cazy_unique_HMM_expression <- na.omit(cazy_unique_HMM_expression)

# Get number of HMMs for plotting (calculating height)
cazy_unique_HMM_expression %>% 
  group_by(family) %>% 
  summarise(number = length(unique(HMM))) -> number_family
number_family <- as.data.frame(number_family)

###################
## GH
###################

pdf(file = paste0(mydir, "/Figures/GH.pdf"), width = 5, height = (number_family[number_family$family == "GH",2])/2)
cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "GH") %>% 
  select(!family) %>% 
  count(HMM, "GH", name = "all") %>% 
  select(!'"GH"') -> GH_1

cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "GH" & up == TRUE) %>% 
  select(!family) %>% 
  count(HMM, "GH", name = "up") %>% 
  select(!'"GH"') %>% 
  right_join(GH_1, by = c("HMM")) %>%                 
  arrange(-all) %>% 
  ggplot(aes(y = fct_reorder(HMM, all), x = all)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_bar(stat = "identity", width = 0.8, aes(x = up), fill = "red") +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(position = "top", name = "Number of genes") +
  ggtitle("GH") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        plot.title = element_text(size = 48, hjust = 0.5))
dev.off()                         

###################
## AA
###################

pdf(file = paste0(mydir, "/Figures/AA.pdf"), width = 5, height = (number_family[number_family$family == "AA",2])/2)
cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "AA") %>% 
  select(!family) %>% 
  count(HMM, "AA", name = "all") %>% 
  select(!'"AA"') -> AA_1

cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "AA" & up == TRUE) %>% 
  select(!family) %>% 
  count(HMM, "AA", name = "up") %>% 
  select(!'"AA"') %>% 
  right_join(AA_1, by = c("HMM")) %>%                 
  arrange(-all) %>% 
  ggplot(aes(y = fct_reorder(HMM, all), x = all)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_bar(stat = "identity", width = 0.8, aes(x = up), fill = "red") +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(position = "top", name = "Number of genes") +
  ggtitle("AA") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        plot.title = element_text(size = 48, hjust = 0.5))
dev.off()

###################
## PL
###################

pdf(file = paste0(mydir, "/Figures/PL.pdf"), width = 5, height = 2)
cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "PL") %>% 
  select(!family) %>% 
  count(HMM, "PL", name = "all") %>% 
  select(!'"PL"') -> PL_1

cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "PL" & up == TRUE) %>% 
  select(!family) %>% 
  count(HMM, "PL", name = "up") %>% 
  select(!'"PL"') %>% 
  right_join(PL_1, by = c("HMM")) %>%                 
  arrange(-all) %>% 
  ggplot(aes(y = fct_reorder(HMM, all), x = all)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_bar(stat = "identity", width = 0.8, aes(x = up), fill = "red") +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(position = "top", name = "Number of genes") +
  ggtitle("PL") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        plot.title = element_text(size = 48, hjust = 0.5))
dev.off()

###################
## CE
###################

pdf(file = paste0(mydir, "/Figures/CE.pdf"), width = 5, height = (number_family[number_family$family == "CE",2])/2)
cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "CE") %>% 
  select(!family) %>% 
  count(HMM, "CE", name = "all") %>% 
  select(!'"CE"') -> CE_1

cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "CE" & up == TRUE) %>% 
  select(!family) %>% 
  count(HMM, "CE", name = "up") %>% 
  select(!'"CE"') %>% 
  right_join(CE_1, by = c("HMM")) %>%                 
  arrange(-all) %>% 
  ggplot(aes(y = fct_reorder(HMM, all), x = all)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_bar(stat = "identity", width = 0.8, aes(x = up), fill = "red") +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(position = "top", name = "Number of genes") +
  ggtitle("CE") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        plot.title = element_text(size = 48, hjust = 0.5))
dev.off()


###################
## GT
###################

pdf(file = paste0(mydir, "/Figures/GT.pdf"), width = 5, height = (number_family[number_family$family == "GT",2])/2)
cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "GT") %>% 
  select(!family) %>% 
  count(HMM, "GT", name = "all") %>% 
  select(!'"GT"') -> GT_1

cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "GT" & up == TRUE) %>% 
  select(!family) %>% 
  count(HMM, "GT", name = "up") %>% 
  select(!'"GT"') %>% 
  right_join(GT_1, by = c("HMM")) %>%                 
  arrange(-all) %>% 
  ggplot(aes(y = fct_reorder(HMM, all), x = all)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_bar(stat = "identity", width = 0.8, aes(x = up), fill = "red") +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(position = "top", name = "Number of genes") +
  ggtitle("GT") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        plot.title = element_text(size = 48, hjust = 0.5))
dev.off()

###################
## CBM
###################

pdf(file = paste0(mydir, "/Figures/CBM.pdf"), width = 5, height = (number_family[number_family$family == "CBM",2])/2)
cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "CBM") %>% 
  select(!family) %>% 
  count(HMM, "CBM", name = "all") %>% 
  select(!'"CBM"') -> CBM_1

cazy_unique_HMM_expression %>%
  ungroup() %>% 
  select(!gene) %>% 
  filter(family == "CBM" & up == TRUE) %>% 
  select(!family) %>% 
  count(HMM, "CBM", name = "up") %>% 
  select(!'"CBM"') %>% 
  right_join(CBM_1, by = c("HMM")) %>%                 
  arrange(-all) %>% 
  ggplot(aes(y = fct_reorder(HMM, all), x = all)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_bar(stat = "identity", width = 0.8, aes(x = up), fill = "red") +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(position = "top", name = "Number of genes") +
  ggtitle("CBM") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        plot.title = element_text(size = 48, hjust = 0.5))
dev.off()


###################
## Dotplot 
###################

cazy_unique_HMM_expression %>% 
  ungroup() %>% 
  left_join(gene_expression, by = "gene") %>% 
  group_by(gene) %>% 
  #mutate(attacking_TPM = mean(c(V1S7_TPM, V1S8_TPM, V1S9_TPM))) %>% 
  ungroup() %>% 
  select(gene, HMM, family, up, down, baseMean) -> dotplot_data

dotplot_data %>% 
  mutate(status = case_when(dotplot_data$up == "TRUE" ~ "up",
                            dotplot_data$down == "TRUE" ~ "down",
                            dotplot_data$up == "FALSE" & dotplot_data$down == "FALSE" ~ "none")) -> dotplot_data

cols <- c("up" = "#DD3C51", "down" = "#5FA1F7", "none" = "#595959")

###################
## CBM
###################

filter(dotplot_data, family == "CBM") %>% 
  group_by(HMM) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) -> CBM_families_ranked 

CBM_families <- CBM_families_ranked$HMM
CBM_plots <- list()
max_genes <- c()

for (i in CBM_families) {
filter(dotplot_data, HMM == i) %>% 
  mutate(new = fct_reorder(gene, baseMean, .desc = TRUE)) %>% 
  ggplot(aes(x = new, y = log2(baseMean), fill = status)) +
    geom_bar(stat = "identity", aes(x = new, y = 19.6359), fill = "grey90") +
    geom_bar(stat = "identity") +
    geom_bar(stat = "identity", aes(x = new, y = 19.6359), fill = NA, color = "black") +
    scale_fill_manual(values = cols, guide = NULL) +
    scale_y_continuous(limits = c(0, 19.6359), 
                       name = NULL, 
                       breaks = NULL,
                       oob = rescale_none) +
    scale_x_discrete(name = NULL, guide = NULL) +
    theme_void() +
    theme(panel.grid = element_blank()) -> g
  CBM_plots[[i]] <- ggplotGrob(g)
  filter(dotplot_data, HMM == i) %>% 
    nrow() -> max_genes[i]
}

# Write script to create plot
param <- c()
for (i in 1:length(CBM_plots)) {
  a <- print(paste0("draw_plot(CBM_plots[[", i, "]], x = x, y = rel_position[",i,"], width = y*(max_genes[",i,"]/max(max_genes)), height = h) +"))
  param[i] <- a
}

sink(file = paste0(mydir, "/Figures/Figures_cazy/CBM_plot.R"))
cat(paste0("rel_position <- seq(0.9, 0.15, -0.05)", "\n"),
    paste0("pdf(file = paste0(mydir, \"/Figures/CAZY_overview_CBM.pdf\"), width = 15, height = 15)", "\n"),
    paste0("h = 0.035", "\n"),
    paste0("x = 0.1", "\n"),
    paste0("y = 0.8", "\n"),
    paste0("ggdraw() +", "\n"), 
    paste0(param, "\n"), 
    paste0("draw_plot_label(label = CBM_families, x = c(0), y = c(rel_position + 0.025), size = 12)", "\n"),
    paste0("dev.off()"))
sink()

###################
## Cleaving enzymes
###################

filter(dotplot_data, family == "GH") %>% 
  group_by(HMM) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) -> GH_families_ranked 

filter(dotplot_data, HMM == "AA11") %>% 
  group_by(HMM) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) -> AA_families_ranked 

filter(dotplot_data, family == "PL") %>% 
  group_by(HMM) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) -> PL_families_ranked 

cleaving_families_ranked <- rbind(GH_families_ranked,
                                  PL_families_ranked,
                                  AA_families_ranked)

cleaving_families <- cleaving_families_ranked$HMM
cleaving_plots <- list()
max_genes <- c()

for (i in cleaving_families) {
    filter(dotplot_data, HMM == i) %>% 
      mutate(new = fct_reorder(gene, baseMean, .desc = TRUE)) %>% 
      ggplot(aes(x = new, y = log2(baseMean), fill = status)) +
        geom_bar(stat = "identity", aes(x = new, y = 19.6359), fill = "grey90") +
        geom_bar(stat = "identity") +
        geom_bar(stat = "identity", aes(x = new, y = 19.6359), fill = NA, color = "black") +
        scale_fill_manual(values = cols, guide = NULL) +
        scale_y_continuous(limits = c(0, 19.6359), 
                          name = NULL, 
                          breaks = NULL,
                          oob = rescale_none) +
        scale_x_discrete(name = NULL, guide = NULL) +
        theme_void() +
        theme(panel.grid = element_blank()) -> g
    cleaving_plots[[i]] <- ggplotGrob(g)
    filter(dotplot_data, HMM == i) %>% 
      nrow() -> max_genes[i]
}

# Write script to create plot
param <- c()
for (i in 1:length(cleaving_plots)) {
  a <- print(paste0("draw_plot(cleaving_plots[[", i, "]], x = x, y = rel_position[",i,"], width = y*(max_genes[",i,"]/max(max_genes)), height = h) +"))
  param[i] <- a
}

sink(file = paste0(mydir, "/Figures/Figures_cazy/cleaving_plot.R"))
cat(paste0("rel_position <- seq(0.95, 0.1, -0.0131)", "\n"),
    paste0("pdf(file = paste0(mydir, \"/Figures/CAZY_overview_cleaving.pdf\"), width = 15, height = 30)", "\n"),
    paste0("h = 0.01", "\n"),
    paste0("x = 0.07", "\n"),
    paste0("y = 0.2", "\n"),
    paste0("ggdraw() +", "\n"), 
    paste0(param, "\n"), 
    paste0("draw_plot_label(label = cleaving_families, x = c(0), y = c(rel_position + 0.008), size = 8)", "\n"),
    paste0("dev.off()"))
sink()

###################
## Figure ranked TPM, CBM
###################


cazy_unique_HMM_expression %>% 
  ungroup() %>% 
  left_join(gene_expression, by = "gene") %>% 
  group_by(gene) %>% 
  mutate(attacking_TPM = mean(c(V1S7_TPM, V1S8_TPM, V1S9_TPM))) %>% 
  ungroup() %>% 
  select(gene, HMM, family, up, down, attacking_TPM) -> dotplot_data

dotplot_data %>% 
  mutate(status = case_when(dotplot_data$up == "TRUE" ~ "up",
                            dotplot_data$down == "TRUE" ~ "down",
                            dotplot_data$up == "FALSE" & dotplot_data$down == "FALSE" ~ "none")) -> dotplot_data

#cols <- c("up" = "#ca0020", "down" = "#0571b0", "none" = "#595959")

dotplot_data %>% 
  filter(family == "CBM") %>% 
  group_by(gene) %>% 
  distinct(HMM, .keep_all = TRUE) %>% 
  mutate(HMM = paste0(HMM, collapse = "+")) %>%
  distinct(HMM, .keep_all = TRUE) %>% 
  ungroup()  %>% 
  arrange(-attacking_TPM) -> attacking_TPM_df

# Plot the top 50 in decreasing order

#pdf(file = paste0(mydir, "/Figures/CBM_TPM.pdf"), height = 12, width = 3.5)
attacking_TPM_df[c(1:50),] %>% 
  mutate(new = fct_reorder(gene, attacking_TPM)) %>% 
  ggplot(aes(y = new, x = attacking_TPM, fill = status)) +
  geom_col(width = 0.8) +
  scale_y_discrete(name = NULL, labels = rev(attacking_TPM_df$HMM[c(1:50)])) +
  scale_x_continuous(name = "TPM") +
  scale_fill_manual(values = cols, guide = NULL) +
  theme_bw() -> a
#dev.off()

###################
## Figure ranked TPM, cleaving
###################


cazy_unique_HMM_expression %>% 
  ungroup() %>% 
  left_join(gene_expression, by = "gene") %>% 
  group_by(gene) %>% 
  mutate(attacking_TPM = mean(c(V1S7_TPM, V1S8_TPM, V1S9_TPM))) %>% 
  ungroup() %>% 
  select(gene, HMM, family, up, down, attacking_TPM) -> dotplot_data

dotplot_data %>% 
  mutate(status = case_when(dotplot_data$up == "TRUE" ~ "up",
                            dotplot_data$down == "TRUE" ~ "down",
                            dotplot_data$up == "FALSE" & dotplot_data$down == "FALSE" ~ "none")) -> dotplot_data

#cols <- c("up" = "#d7191c", "down" = "#2c7bb6", "none" = "#595959")

filter(dotplot_data, family == "GH") %>% 
  group_by(HMM) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) -> GH_families_ranked 

filter(dotplot_data, HMM == "AA11") %>% 
  group_by(HMM) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) -> AA_families_ranked 

filter(dotplot_data, family == "PL") %>% 
  group_by(HMM) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) -> PL_families_ranked 

cleaving_families_ranked <- rbind(GH_families_ranked,
                                  PL_families_ranked,
                                  AA_families_ranked)

cleaving_families <- cleaving_families_ranked$HMM


dotplot_data %>% 
  filter(HMM %in% cleaving_families) %>% 
  group_by(gene) %>% 
  distinct(HMM, .keep_all = TRUE) %>% 
  mutate(HMM = paste0(HMM, collapse = "+")) %>%
  distinct(HMM, .keep_all = TRUE) %>% 
  ungroup()  %>% 
  arrange(-attacking_TPM) -> attacking_TPM_df

# Plot the top 50 in decreasing order

#pdf(file = paste0(mydir, "/Figures/Lytic_TPM.pdf"), height = 12, width = 3.5)
attacking_TPM_df[c(1:50),] %>% 
  mutate(new = fct_reorder(gene, attacking_TPM)) %>% 
  ggplot(aes(y = new, x = attacking_TPM, fill = status)) +
  geom_col(width = 0.8) +
  scale_y_discrete(name = NULL, labels = rev(attacking_TPM_df$HMM[c(1:50)])) +
  scale_x_continuous(name = "TPM") +
  scale_fill_manual(values = cols, guide = NULL) +
  theme_bw() -> b
#dev.off()

pdf(file = paste0(mydir, "/Figures/TPM.pdf"), width = 8.2, height = 11.6)
grid::grid.newpage()
grid::grid.draw(cbind(ggplotGrob(a), 
                      ggplotGrob(b)))
dev.off()

###################
## Figure 3B
###################

# Rank CAZY-annotated proteins according to TPM in each condition

cazy_expression <- read.delim(file = paste0(mydir, "/Tables/cazy.expression.out"))
cazy_expression$HMMER <- gsub("\\(.*?\\)", "", cazy_expression$HMMER)

seeking <- c("V1S4_TPM", "V1S5_TPM", "V1S6_TPM")
attacking <- c("V1S7_TPM", "V1S8_TPM", "V1S9_TPM")
digesting <- c("V1S1_TPM", "V1S2_TPM", "V1S3_TPM")

cazy_expression %>% 
  select(c(1, 3:11, 18, 19, 20)) %>% 
  pivot_longer(cols = starts_with("V1S"), names_to = "sample", values_to = "TPM") %>% 
  mutate(group = case_when(sample %in% seeking ~ "seeking",
                           sample %in% attacking ~ "attacking",
                           sample %in% digesting ~ "digesting")) %>% 
  group_by(protein, group) %>% 
  mutate(mean_TPM = mean(TPM)) %>% 
  ungroup() %>% 
  mutate(family = case_when(grepl("GH", HMMER) ~ "GH",
                            grepl("GT", HMMER) ~ "GT",
                            grepl("PL", HMMER) ~ "PL",
                            grepl("CE", HMMER) ~ "CE",
                            grepl("AA", HMMER) ~ "AA",
                            grepl("CBM", HMMER) ~ "CBM")) -> cazy_expression_longer

# Seeking

cazy_expression_longer %>% 
  filter(group == "seeking") %>% 
  select(!c(5,6,7)) %>% 
  distinct(protein, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  distinct(HMMER, .keep_all = TRUE) %>% 
  mutate(HMM = paste0(HMMER, collapse = "+")) %>%
  distinct(HMM, .keep_all = TRUE) %>% 
  ungroup()  %>% 
  arrange(-mean_TPM) -> seeking_df

# Plot the top 50 in decreasing order

seeking_df[c(1:50),] %>% 
  mutate(new = fct_reorder(protein, mean_TPM)) %>% 
  ggplot(aes(y = new, x = mean_TPM, fill = family)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(name = NULL, labels = rev(seeking_df$HMM[c(1:50)])) +
    scale_x_continuous(name = "TPM") +
    scale_fill_ochre() +
    theme_bw()

# Attacking

cazy_expression_longer %>% 
  filter(group == "attacking") %>% 
  select(!c(5,6,7)) %>% 
  distinct(protein, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  distinct(HMMER, .keep_all = TRUE) %>% 
  mutate(HMM = paste0(HMMER, collapse = "+")) %>%
  distinct(HMM, .keep_all = TRUE) %>% 
  ungroup()  %>% 
  arrange(-mean_TPM) -> attacking_df

# Plot the top 50 in decreasing order

attacking_df[c(1:50),] %>% 
  mutate(new = fct_reorder(protein, mean_TPM)) %>% 
  ggplot(aes(y = new, x = mean_TPM, fill = family)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(name = NULL, labels = rev(attacking_df$HMM[c(1:50)])) +
  scale_x_continuous(name = "TPM") +
  scale_fill_ochre() +
  theme_bw()

# Digesting

cazy_expression_longer %>% 
  filter(group == "digesting") %>% 
  select(!c(5,6,7)) %>% 
  distinct(protein, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  distinct(HMMER, .keep_all = TRUE) %>% 
  mutate(HMM = paste0(HMMER, collapse = "+")) %>%
  distinct(HMM, .keep_all = TRUE) %>% 
  ungroup()  %>% 
  arrange(-mean_TPM) -> digesting_df

# Plot the top 50 in decreasing order

digesting_df[c(1:50),] %>% 
  mutate(new = fct_reorder(protein, mean_TPM)) %>% 
  ggplot(aes(y = new, x = mean_TPM, fill = family)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(name = NULL, labels = rev(digesting_df$HMM[c(1:50)])) +
  scale_x_continuous(name = "TPM") +
  scale_fill_ochre() +
  theme_bw()
  
###################
## Figure MEROPS
###################

merops <- read.delim(file = paste0(mydir, "/Tables/merops.expression.gene.out"))
merops %>% 
  select(c(1, 2, 26, 27, 24)) -> merops
merops$annotation <- paste(merops$gene, merops$CATALYTIC.TYPE, merops$TYPE.ENZYME, merops$id, sep = ", ")
merops <- arrange(merops, -baseMean)
merops[c(1:50),] %>% 
  mutate(new = fct_reorder(annotation, baseMean)) %>% 
  ggplot(aes(y = new, x = baseMean)) +
  geom_bar(stat = "identity", fill = "#786060") +
  scale_y_discrete(name = NULL, labels = rev(merops$annotation[c(1:50)])) +
  scale_x_continuous(name = "baseMean") +
  theme_bw()

' Figure 3A without gene expression information
pdf(file = paste0(mydir, "/Figures/GH.pdf"), width = 5, height = (number_family[number_family$family == "GH",2])/2)
ggplot(filter(cazy_unique_HMM, family == "GH"), aes(y = fct_rev(fct_infreq(HMM)))) +
  geom_bar(width = 0.8) +
  geom_bar()
scale_y_discrete(name = "GH") +
  scale_x_continuous(position = "top", name = "Number of genes") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()

pdf(file = paste0(mydir, "/Figures/PL.pdf"), width = 5, height = (number_family[number_family$family == "PL",2])/2)
ggplot(filter(cazy_unique_HMM, family == "PL"), aes(y = fct_rev(fct_infreq(HMM)))) +
  geom_bar() +
  scale_y_discrete(name = "PL") +
  scale_x_continuous(position = "top", name = "Number of genes") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()

pdf(file = paste0(mydir, "/Figures/CBM.pdf"), width = 5, height = (number_family[number_family$family == "CBM",2])/2)
ggplot(filter(cazy_unique_HMM, family == "CBM"), aes(y = fct_rev(fct_infreq(HMM)))) +
  geom_bar(width = 0.8) +
  scale_y_discrete(name = "CBM") +
  scale_x_continuous(position = "top", name = "Number of genes") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()

pdf(file = paste0(mydir, "/Figures/CE.pdf"), width = 5, height = (number_family[number_family$family == "CE",2])/2)
ggplot(filter(cazy_unique_HMM, family == "CE"), aes(y = fct_rev(fct_infreq(HMM)))) +
  geom_bar() +
  scale_y_discrete(name = "CE") +
  scale_x_continuous(position = "top", name = "Number of genes") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()

pdf(file = paste0(mydir, "/Figures/GT.pdf"), width = 5, height = (number_family[number_family$family == "GT",2])/2)
ggplot(filter(cazy_unique_HMM, family == "GT"), aes(y = fct_rev(fct_infreq(HMM)))) +
  geom_bar() +
  scale_y_discrete(name = "GT") +
  scale_x_continuous(position = "top", name = "Number of genes") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()

pdf(file = paste0(mydir, "/Figures/AA.pdf"), width = 5, height = (number_family[number_family$family == "AA",2])/2)
ggplot(filter(cazy_unique_HMM, family == "AA"), aes(y = fct_rev(fct_infreq(HMM)))) +
  geom_bar() +
  scale_y_discrete(name = "AA") +
  scale_x_continuous(position = "top", name = "Number of genes") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()'
