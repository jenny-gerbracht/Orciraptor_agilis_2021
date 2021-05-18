library(seqinr)
library(tximport)
library(dplyr)
library(stringi)
library(stringr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(Glimma)
library(ComplexHeatmap)
library(eply)
library(RColorBrewer)
library(circlize)
library(rhmmer)
library(dendextend)
library(tidyr)
library(grid)
library(gridExtra)
library(forcats)
library(ggpubr)
library(ggrepel)

source(file = "../../config.txt")

###################
## Create gene trans map
###################

# A gene_trans_maps from Trinity has the order gene_id, transcript_id
# tximport and Lace need the map in the order transcript_id, gene_id

assembly <- read.fasta(file = paste0(mydir, "/Module_3/orciraptor_200_filtered.fasta"))
gene_trans_map <- as.data.frame(matrix(0, 
                                       ncol = 2, 
                                       nrow = length(assembly)))
gene_trans_map$V1 <- attr(assembly, "names")
gene_trans_map$V2 <- str_split(gene_trans_map$V1, "_", simplify = TRUE)[,7]
write.table(gene_trans_map, 
            file = paste0(mydir, "/Module_5/gene_trans_map"), 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

# Trinity style
gene_trans_map_trinity <- as.data.frame(matrix(0, 
                                       ncol = 2, 
                                       nrow = length(assembly)))
gene_trans_map_trinity$V2 <- attr(assembly, "names")
gene_trans_map_trinity$V1 <- str_split(gene_trans_map_trinity$V2, "_", simplify = TRUE)[,7]
write.table(gene_trans_map_trinity, 
            file = paste0(mydir, "/Module_5/gene_trans_map_trinity"), 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE, 
            sep = "\t")


###################
## Import expression data from salmon
###################


samples_j <- read.table(paste0(mydir, "/experiment.txt"), 
                        col.names = c("samples", "condition"))
samples_j <- samples_j[c(4,5,6,7,8,9,1,2,3),]
files_j <- file.path(paste0(mydir, "/Module_5"), 
                     "salmon", 
                     paste0(samples_j$samples, ".salmon_quant"), 
                     "quant.sf")
names(files_j) <- samples_j$samples
tx2gene_j <- read.table(file = paste0(mydir, "/Module_5/gene_trans_map"))
all(file.exists(files_j))
txi_j <- tximport(files_j, 
                  type = "salmon", 
                  tx2gene = tx2gene_j)

txi_i <- tximport(files_j, 
                  type = "salmon",
                  txOut = TRUE)

# Export count and abundance (TPM) tables, and isoform TPM table for GO analysis

write.table(txi_j$counts,
            file = paste0(mydir, "/Module_5/salmon_gene_count_matrix"),
            quote = FALSE,
            sep = "\t")

write.table(txi_j$abundance,
            file = paste0(mydir, "/Module_5/salmon_gene_TPM_matrix"),
            quote = FALSE,
            sep = "\t")

write.table(txi_i$abundance,
            file = paste0(mydir, "/Module_5/salmon_isoform_TPM_matrix"),
            quote = FALSE,
            sep = "\t")

###################
## DESeq2
###################


samples_j$condition <- factor(samples_j$condition)
dds <- DESeqDataSetFromTximport(txi_j,
                                colData = samples_j,
                                design = ~ condition)

# Pre-filtering low count genes: analogous to edgeR pre-filtering rowSums(cpm(counts)>1)>=2, gives same results as when
# using keep <- rowSums(counts(dds)) >= 10 which is described in the pre-filtering section of the DESeq2 vignette

dds <- dds[rowSums(fpm(dds) > 1) >= 2] 

dds <- DESeq(dds)
res_as <- results(dds, contrast = c("condition", "attacking", "seeking"))
res_ds <- results(dds, contrast = c("condition", "digesting", "seeking"))
res_da <- results(dds, contrast = c("condition", "digesting", "attacking"))

head(as.data.frame(res_as))

# Extract gene names of all included genes for background set in GO analysis
write.table(as.data.frame(rownames(dds)),
            file = paste0(mydir, "/Module_5/background.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

###################
## PCA plot
###################

# Extracting transformed count values

vsd <- vst(dds, blind = FALSE)

# Make plot

pcaData <- plotPCA(vsd, 
                   intgroup = c("condition"), 
                   returnData=TRUE)
pcaData$condition <- factor(pcaData$condition, levels = c("seeking", "attacking", "digesting"))   
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file = paste0(mydir, "/Figures/PCA.pdf"), width = 5, height = 7)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=c("#3333CC", "#D60093", "#1F4E79")) +
  theme_minimal() +
  theme(panel.border = element_rect(fill = "NA"))
dev.off()

###################
## Building output tables
###################

# Base table gene expression info
# TPM values
gene_expression <- as.data.frame(txi_j$abundance)
colnames(gene_expression) <- paste0(colnames(gene_expression), "_TPM")
gene_expression <- rownames_to_column(gene_expression, var = "gene")

# Add baseMean
as.data.frame(res_as) %>%
  rownames_to_column(var = "gene") %>% 
  select(gene, baseMean) %>% 
  left_join(gene_expression, by = "gene") -> gene_expression

# Add log2FC and padj values from each comparison
as.data.frame(res_as) %>% 
  setNames(., paste0(colnames(res_as), "_as")) %>% 
  rownames_to_column(var = "gene") %>% 
  select(gene, log2FoldChange_as, padj_as) %>% 
  left_join(gene_expression, by = "gene") -> gene_expression

as.data.frame(res_da) %>% 
  setNames(., paste0(colnames(res_da), "_da")) %>% 
  rownames_to_column(var = "gene") %>% 
  select(gene, log2FoldChange_da, padj_da) %>% 
  left_join(gene_expression, by = "gene") -> gene_expression

as.data.frame(res_ds) %>% 
  setNames(., paste0(colnames(res_ds), "_ds")) %>% 
  rownames_to_column(var = "gene") %>% 
  select(gene, log2FoldChange_ds, padj_ds) %>% 
  left_join(gene_expression, by = "gene") -> gene_expression

gene_expression %>% 
  select(gene, baseMean, V1S4_TPM, V1S5_TPM, V1S6_TPM, V1S7_TPM, V1S8_TPM, V1S9_TPM, V1S1_TPM, V1S2_TPM, V1S3_TPM, log2FoldChange_as, padj_as, log2FoldChange_da, padj_da, log2FoldChange_ds, padj_ds) -> gene_expression

write.table(gene_expression, 
            file = paste0(mydir, "/Tables/gene.expression.out"), 
            quote = FALSE, 
            sep = "\t")

###################
## Clustering
###################

# Cutoff log2FC 1, padj 0.05

as.data.frame(res_as) %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) -> res_as_2
as.data.frame(res_da) %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) -> res_da_2
as.data.frame(res_ds) %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) -> res_ds_2
dge_set_2 <- c(rownames(res_as_2), rownames(res_da_2), rownames(res_ds_2))
dge_set_2 <- unique(dge_set_2)
#Generate z-score matrix for heatmaps
vsd_matrix <- assay(vsd)
vsd_matrix_dge_2 <- vsd_matrix[dge_set_2,]
heat_2 <- t(scale(t(vsd_matrix_dge_2)))

# Clusters rows by Pearson correlation as distance method
hc <- hclust(as.dist(1 - cor(t(as.matrix(vsd_matrix_dge_2)))))
plot(as.dendrogram(hc), leaflab = "none")
my_gene_partition_assignments2 <- cutree(hc, h=80/100*max(hc$height))

# Visualise dendrogram with clusters

clust.cutree <- dendextend:::cutree(as.dendrogram(hc),h=80/100*max(hc$height), order_clusters_as_data = FALSE)

idx <- order(names(clust.cutree))
clust.cutree <- clust.cutree[idx]
df.merge <- merge(my_gene_partition_assignments2,clust.cutree,by='row.names')
df.merge.sorted <- df.merge[order(df.merge$y),]
lbls<-unique(df.merge.sorted$x)
dend1 <- color_branches(as.dendrogram(hc),h=80/100*max(hc$height), groupLabels = lbls)
plot(dend1)
abline(h=80/100*max(hc$height), lty = 2, col="grey")

# Make factor_labeling.txt for goseq
data.frame(my_gene_partition_assignments2) %>% 
  rownames_to_column(var = "genes") -> factor_labeling
factor_labeling$my_gene_partition_assignments2 <- factor_labeling$my_gene_partition_assignments
write.table(factor_labeling, 
            file = paste0(mydir, "/Module_5/GO_analysis/factor_labeling.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t")

# Make list of clusters
clusterlist <- list()
for (i in c(1:max(my_gene_partition_assignments2))) {
  cluster <- heat_2[(my_gene_partition_assignments2 == i),]
  clusterlist[[i]] <- cluster
}

# Plot clusters as boxplots
p <- list()
splan <- 3 - 1L

for (i in 1:max(my_gene_partition_assignments2)) {
  box <- rownames_to_column(as.data.frame(clusterlist[[i]]), var = "gene")
  box_longer <- pivot_longer(box, cols = !gene, names_to = "sample")
  box_longer %>% 
    filter(sample == "V1S4" | sample == "V1S5" | sample == "V1S6") %>% 
    mutate(cat = "Seeking") -> box_longer1
  box_longer %>% 
    filter(sample == "V1S7" | sample == "V1S8" | sample == "V1S9") %>% 
    mutate(cat = "Attacking") -> box_longer2
  box_longer %>% 
    filter(sample == "V1S1" | sample == "V1S2" | sample == "V1S3") %>% 
    mutate(cat = "Digesting") -> box_longer3
  comb <- rbind(box_longer1, box_longer2, box_longer3)
  comb$cat <-factor(comb$cat, levels = c("Seeking", "Attacking", "Digesting"))
  g <- ggplot(comb, aes(x = cat, y = value)) +
    geom_boxplot(outlier.size = 0,
                 outlier.shape = NA,
                 alpha = 0.5, 
                 color = "#252A52",
                 fill = "#252A52") +
    stat_smooth(aes(x = cat, y = value, group = "#252A52"), color = "#252A52",
                se = TRUE,
                method = "lm", formula = y~poly(x, splan)) +
    ggtitle(paste("cluster", i, "|", "genes:", nrow(clusterlist[[i]]))) +
    scale_y_continuous(limits = c(-2, 2)) +
    ylab("z-score") +
    xlab(NULL) +
    theme(legend.position="none", panel.background = element_rect(colour = "darkgrey", size=1))
  p[[i]] <- ggplotGrob(g)
}
rm(box, box_longer, box_longer1, box_longer2, box_longer3, comb)
pdf(file = paste0(mydir, "/Figures/cluster_boxplots_log2FC1.pdf"), width = 15, height = 3)
grid.arrange(grobs = p, ncol=5)
dev.off()



###################
## Functional annotation
###################
# MEROPS
# Output table for all peptide sequences

merops_out <- read.delim(paste0(mydir, "/Module_4/merops/merops.out"), 
                         header = FALSE,
                         col.names = c("qseqid", "sseqid", "stitle", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
merops_out$gene <- str_split(merops_out$qseqid, "_", simplify = TRUE)[,1]
merops_out %>% 
  select(gene, qseqid, sseqid, stitle, pident, length, evalue) %>% 
  left_join(gene_expression, by = "gene") -> merops_expression
merops_expression$isoform <- str_split(merops_expression$qseqid, "\\.", simplify = TRUE)[,1]
merops_expression <- rename(merops_expression, "qseqid" = "protein")
merops_expression %>% 
  select(gene, baseMean, V1S4_TPM, V1S5_TPM, V1S6_TPM, V1S7_TPM, V1S8_TPM, V1S9_TPM, V1S1_TPM, V1S2_TPM, V1S3_TPM, log2FoldChange_as, padj_as, log2FoldChange_da, padj_da, log2FoldChange_ds, padj_ds,
         isoform, protein, stitle, pident, length, evalue) -> merops_expression
write.table(merops_expression, 
            file = paste0(mydir, "/Tables/merops.expression.out"), 
            quote = FALSE, 
            sep = "\t")

# Keep entry per gene with lowest evalue, extract merops family assignment, match by family

merops_out <- read.delim(paste0(mydir, "/Module_4/merops/merops.out"), 
                         header = FALSE,
                         col.names = c("qseqid", "sseqid", "stitle", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
merops_out$gene <- str_split(merops_out$qseqid, "_", simplify = TRUE)[,1]
merops_out %>% 
  select(gene, qseqid, sseqid, stitle, pident, length, evalue) %>% 
  left_join(gene_expression, by = "gene") -> merops_expression
merops_expression$isoform <- str_split(merops_expression$qseqid, "\\.", simplify = TRUE)[,1]
merops_expression <- rename(merops_expression, "qseqid" = "protein")
merops_expression %>% 
  select(gene, baseMean, V1S4_TPM, V1S5_TPM, V1S6_TPM, V1S7_TPM, V1S8_TPM, V1S9_TPM, V1S1_TPM, V1S2_TPM, V1S3_TPM, log2FoldChange_as, padj_as, log2FoldChange_da, padj_da, log2FoldChange_ds, padj_ds,
         isoform, protein, stitle, pident, length, evalue) %>% 
  group_by(gene) %>% 
  mutate(min_evalue = min(evalue)) %>% 
  filter(evalue == min_evalue) %>% 
  group_by(gene) %>% 
  distinct(stitle, .keep_all = TRUE) ->  merops_expression_gene
merops_expression_gene$stitle <- str_replace(merops_expression_gene$stitle, "\\[misleading\\]", "")
merops_expression_gene$id <- str_extract_all(merops_expression_gene$stitle, "(?<=\\[).+?(?=\\])", simplify = TRUE)
merops_expression_gene$FAMILY <- str_split(merops_expression_gene$id, "\\.", simplify = TRUE)[,1]
merops_families <- read.table(paste0(mydir, "/Tables/MEROPS_families.txt"), 
                              sep = "\t", 
                              header = TRUE)
merops_expression_gene %>% 
  left_join(merops_families, by = "FAMILY")  %>% 
  select(!min_evalue) -> merops_expression_gene


write.table(merops_expression_gene, 
            file = paste0(mydir, "/Tables/merops.expression.gene.out"), 
            quote = FALSE, 
            sep = "\t")

# CAZY

cazy_out <- read.delim(paste0(mydir, "/Module_4/cazy_05/overview.txt"))
cazy_out$gene <- str_split(cazy_out$Gene.ID, "_", simplify = TRUE)[,1]
cazy_out %>% 
  select(gene, Gene.ID, HMMER) %>% 
  left_join(gene_expression, by = "gene") -> cazy_out
cazy_out$isoform <- str_split(cazy_out$Gene.ID, "\\.", simplify = TRUE)[,1]
cazy_out <- rename(cazy_out, "Gene.ID" = "protein")
cazy_out %>% 
  select(gene, baseMean, V1S4_TPM, V1S5_TPM, V1S6_TPM, V1S7_TPM, V1S8_TPM, V1S9_TPM, V1S1_TPM, V1S2_TPM, V1S3_TPM, log2FoldChange_as, padj_as, log2FoldChange_da, padj_da, log2FoldChange_ds, padj_ds,
         isoform, protein, HMMER) -> cazy_out
write.table(cazy_out, 
            file = paste0(mydir, "/Tables/cazy.expression.out"), 
            quote =   FALSE, 
            sep = "\t")

# eggnog
# diamond

eggnog_diamond <- read.delim(paste0(mydir, "/Module_4/eggnog/orciraptor_diamond.emapper.annotations"), 
                             skip = 4,
                             quote = "")
eggnog_diamond <- eggnog_diamond[-c((nrow(eggnog_diamond)-2):(nrow(eggnog_diamond))),]
eggnog_diamond$source <- c("diamond")

# HMM
eggnog_hmm <- read.delim(paste0(mydir, "/Module_4/eggnog/orciraptor_hmm.emapper.annotations"), 
                         skip = 4,
                         quote = "")
eggnog_hmm <- eggnog_hmm[-c((nrow(eggnog_hmm)-2):(nrow(eggnog_hmm))),]
eggnog_hmm$source <- c("hmm")

# Keep all annotations from Diamond, add annotations for peptide sequences that were annotated only in HMM mode
eggnog_combined <- rbind(eggnog_diamond, (eggnog_hmm[!eggnog_hmm$X.query_name %in% eggnog_diamond$X.query_name,]))

#Extract K numbers for KEGG pathway mapping
eggnog_combined %>% 
  select(X.query_name, KEGG_ko) -> eggnog_kegg
eggnog_kegg$KEGG_ko <- gsub(eggnog_kegg$KEGG_ko, pattern = "ko:", replacement = "")  
eggnog_kegg <- separate_rows(eggnog_kegg, KEGG_ko, sep = ",")
write.table(eggnog_kegg, 
            file = paste0(mydir, "/Tables/eggnog_kegg"), 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

###################
## Glimma plots
###################

# attacking vs starving

merops_results_as <- res_as[rownames(res_as) %in% merops_expression_gene$gene,]
status_as <- as.numeric(merops_results_as$padj < 0.05 & abs(merops_results_as$log2FoldChange) > 1)
anno_as <- data.frame(gene = rownames(merops_results_as))
anno_as %>% 
  left_join(select(merops_expression_gene, gene, CATALYTIC.TYPE, TYPE.ENZYME, id), by = "gene") -> anno_as
colnames(anno_as) <- c("GeneID", "CATALYTIC.TYPE", "TYPE.ENZYME", "id") 
glMDPlot(merops_results_as, 
         status = status_as, 
         counts = log(counts(dds, normalized = TRUE)[rownames(dds) %in% merops_expression_gene$gene,] + 0.5),
         groups = dds$condition,
         transform = FALSE,
         samples = colnames(dds),
         anno = anno_as,
         path = (paste0(mydir, "/Tables/Glimma")),
         folder = "MEROPS",
         launch = FALSE)


###################
## Volcano plot
###################

cazy_out %>% 
  select(gene, HMMER) %>% 
  mutate(HMMER = gsub("(\\(.*?\\))", "", cazy_out$HMMER)) %>% 
  group_by(gene) %>% 
  distinct(HMMER, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  summarise(HMMER = paste0(HMMER, collapse = "+")) -> genes_cazy


genes_cazy$HMMER <- str_split(genes_cazy$HMMER, "\\+", simplify = TRUE)[,1]
genes_cazy$HMMER <- paste0(genes_cazy$HMMER, " (", genes_cazy$gene, ")")

as.data.frame(res_as) %>% 
  rownames_to_column(var = "gene") -> dge_as
dge_as$cazy <- ifelse(dge_as$gene %in% genes_cazy$gene, "TRUE", "FALSE")
dge_as %>% 
  left_join(genes_cazy) -> dge_as
dge_as$HMMER[dge_as$cazy == FALSE | 
               dge_as$log2FoldChange < 1.5 |
               dge_as$padj > 0.05] <- ""


ggplot(data = dge_as, 
       aes(x = log2FoldChange, 
           y = -log10(pvalue),
           col = cazy, 
           alpha = cazy,
           label = HMMER)) +
  geom_point() +
  geom_point(data = filter(dge_as, cazy == TRUE)) +
  geom_text_repel(color = "black",
                  max.overlaps = Inf) +
  scale_color_manual(values = c("grey", "royalblue")) +
  scale_alpha_manual(values = c(0.2, 0.8)) +
  theme_minimal() +
  geom_vline(xintercept=c(-1.5, 1.5), col = "black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col = "black", linetype = "dashed")

###################
## MA plot
###################

# CAZY de genes with log2(dge_as$baseMean)+1) expression > 12, lfc > 1, padj < 0.05

cazy_out %>% 
  select(gene, HMMER) %>% 
  mutate(HMMER = gsub("(\\(.*?\\))", "", cazy_out$HMMER)) %>% 
  group_by(gene) %>% 
  distinct(HMMER, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  summarise(HMMER = paste0(HMMER, collapse = "+")) -> genes_cazy


genes_cazy$HMMER <- str_split(genes_cazy$HMMER, "\\+", simplify = TRUE)[,1]
genes_cazy$HMMER <- paste0(genes_cazy$HMMER, " (", genes_cazy$gene, ")")

as.data.frame(res_as) %>% 
  rownames_to_column(var = "gene") -> dge_as
dge_as$cazy <- ifelse(dge_as$gene %in% genes_cazy$gene, "TRUE", "FALSE")
dge_as %>% 
  left_join(genes_cazy) -> dge_as
dge_as$HMMER[dge_as$cazy == FALSE | 
               dge_as$log2FoldChange < 1 |
               dge_as$padj > 0.05 |
               (log2(dge_as$baseMean)+1) < 12] <- ""

ggplot(dge_as, aes(x = (log2(baseMean)+1), y = log2FoldChange, label = HMMER)) +
  geom_point(alpha = 0.2, color = "grey") +
  geom_point(data = filter(dge_as, cazy == TRUE), color = "royalblue") +
  scale_x_continuous(breaks=seq(0, max((log2(dge_as$baseMean)+1)), 2)) + 
  labs(x = "Log2 mean expression", y = "Log2 fold change") +
  geom_hline(yintercept = c(0,  1), linetype = c(1, 2),
             color = c("black", "darkgrey")) +
  geom_hline(yintercept = c(0,  2), linetype = c(1, 2),
             color = c("black",  "darkgrey")) +
  geom_text_repel(color = "black",
                  max.overlaps = Inf,
                  force = 20,
                  segment.color = "darkgrey") 

# CAZY genes with log2(dge_as$baseMean)+1) expression > 16
cazy_out %>% 
  select(gene, HMMER) %>% 
  mutate(HMMER = gsub("(\\(.*?\\))", "", cazy_out$HMMER)) %>% 
  group_by(gene) %>% 
  distinct(HMMER, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  summarise(HMMER = paste0(HMMER, collapse = "+")) -> genes_cazy


genes_cazy$HMMER <- str_split(genes_cazy$HMMER, "\\+", simplify = TRUE)[,1]
genes_cazy$HMMER <- paste0(genes_cazy$HMMER, " (", genes_cazy$gene, ")")

as.data.frame(res_as) %>% 
  rownames_to_column(var = "gene") -> dge_as
dge_as$cazy <- ifelse(dge_as$gene %in% genes_cazy$gene, "TRUE", "FALSE")
dge_as %>% 
  left_join(genes_cazy) -> dge_as
dge_as$HMMER[dge_as$cazy == FALSE | 
               (log2(dge_as$baseMean)+1) < 16] <- ""

ggplot(dge_as, aes(x = (log2(baseMean)+1), y = log2FoldChange, label = HMMER)) +
  geom_point(alpha = 0.2, color = "grey") +
  geom_point(data = filter(dge_as, cazy == TRUE), color = "royalblue") +
  scale_x_continuous(breaks=seq(0, max((log2(dge_as$baseMean)+1)), 2)) + 
  labs(x = "Log2 mean expression", y = "Log2 fold change") +
  geom_hline(yintercept = c(0,  1), linetype = c(1, 2),
             color = c("black", "darkgrey")) +
  geom_hline(yintercept = c(0,  2), linetype = c(1, 2),
             color = c("black",  "darkgrey")) +
  geom_text_repel(color = "black",
                  max.overlaps = Inf,
                  force = 15,
                  segment.color = "darkgrey")  
