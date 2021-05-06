library(seqinr)
library(tximport)
library(dplyr)
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

###################
## Create gene trans map
###################

# A gene_trans_maps from Trinity has the order gene_id, transcript_id
# tximport and Lace need the map in the order transcript_id, gene_id

setwd("../../Module_3")
assembly <- read.fasta(file = "orciraptor_200_filtered.fasta")
gene_trans_map <- as.data.frame(matrix(0, ncol = 2, nrow = length(assembly)))
gene_trans_map$V1 <- attr(assembly, "names")
gene_trans_map$V2 <- str_split(gene_trans_map$V1, "_", simplify = TRUE)[,7]
setwd("../Module_5")
write.table(gene_trans_map, 
            file = "gene_trans_map", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

###################
## Import expression data from salmon
###################

dir_j <- getwd()
list.files(dir_j)
samples_j <- read.table(file.path(dir_j, "experiment.txt"), col.names = c("samples", "condition"))
files_j <- file.path(dir_j, "salmon", paste0(samples_j$samples, ".salmon_quant"), "quant.sf")
names(files_j) <- samples_j$samples
tx2gene_j <- read.table("gene_trans_map")
all(file.exists(files_j))
txi_j <- tximport(files_j, 
                  type = "salmon", 
                  tx2gene = tx2gene_j)

# Export count and abundance (TPM) tables

write.table(txi_j$counts,
            file = "salmon_gene_count_matrix",
            quote = FALSE,
            sep = "\t")

write.table(txi_j$abundance,
            file = "salmon_gene_TPM_matrix",
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
setwd("../Figures")
pdf(file="PCA.pdf", width = 5, height = 7)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=c("#3333CC", "#D60093", "#1F4E79")) +
  theme_minimal() +
  theme(panel.border = element_rect(fill = "NA"))
dev.off()

write.table(res_as, file = "gen_names.txt")

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

# MEROPS

setwd("../Module_4/merops")
merops_out <- read.delim("merops.out",
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
setwd("../../Tables/")
write.table(merops_expression, file = "merops.expression.out", quote = FALSE, sep = "\t", )

# CAZY

setwd("../Module_4/cazy_05")
cazy_out <- read.delim("overview.txt")
cazy_out$gene <- str_split(cazy_out$Gene.ID, "_", simplify = TRUE)[,1]
cazy_out %>% 
  select(gene, Gene.ID, HMMER) %>% 
  left_join(gene_expression, by = "gene") -> cazy_out
cazy_out$isoform <- str_split(cazy_out$Gene.ID, "\\.", simplify = TRUE)[,1]
cazy_out <- rename(cazy_out, "Gene.ID" = "protein")
cazy_out %>% 
  select(gene, baseMean, V1S4_TPM, V1S5_TPM, V1S6_TPM, V1S7_TPM, V1S8_TPM, V1S9_TPM, V1S1_TPM, V1S2_TPM, V1S3_TPM, log2FoldChange_as, padj_as, log2FoldChange_da, padj_da, log2FoldChange_ds, padj_ds,
         isoform, protein, HMMER) -> cazy_out
setwd("../../Tables/")
write.table(cazy_out, file = "cazy.expression.out", quote =   FALSE, sep = "\t")

# eggnog
# diamond

setwd("../Module_4/eggnog")
eggnog_diamond <- read.delim("orciraptor_diamond.emapper.annotations", skip = 4)
eggnog_diamond <- eggnog_diamond[-c((nrow(eggnog_diamond)-2):(nrow(eggnog_diamond))),]
