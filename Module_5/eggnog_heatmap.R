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

setwd("~/transfer/Jenny/210104_rnaspades_DGE")

#Prepare tximport
dir_j <- "~/transfer/Jenny/210104_rnaspades_DGE"
list.files(dir_j)
samples_j <- read.table("quant_files.txt")
files_j <- file.path(dir_j, "salmon", samples_j$V1, "quant.sf")
names(files_j) <- read.table("experiment.txt")$V1
tx2gene_j <- read.table("gene_trans_map")[c("V2","V1")]
all(file.exists(files_j))
#Gene level
txi_j <- tximport(files_j, type = "salmon", tx2gene = tx2gene_j)

#Start DESeq2
samples_j2 <- read.delim("experiment.txt", row.names = NULL, sep = "\t",header = FALSE, col.names = c("samples", "condition"))
samples_j2$condition <- sub("control", "seeking", samples_j2$condition)
colnames(samples_j2)
samples_j2$condition <- factor(samples_j2$condition)

dds <- DESeqDataSetFromTximport(txi_j,
                                colData = samples_j2,
                                design = ~ condition)

dds <- dds[rowSums(fpm(dds)>1)>=2] 

dds <- DESeq(dds)
res_fs <- results(dds, contrast = c("condition", "feeding", "seeking"))
res_ds <- results(dds, contrast = c("condition", "digesting", "seeking"))
res_df <- results(dds, contrast = c("condition", "digesting", "feeding"))

#make maplot: cellulases + pectin hydrolases

#Import CAZY annotation
setwd("X:/210129_DGE_logFC1_DESeq2/Annotation")
cazy <- read_tblout("CAZY_noEval.out")
cazy %>% 
  filter(best_domain_evalue < 1E-4) %>% 
  dplyr::select(domain_name, query_name, sequence_evalue) -> cazy_05

#Add gene columns
cazy_05$gene <- sub("_i[0-9][0-9]?[0-9]?.p.", "", cazy_05$query_name)
cazy_05 %>% 
  mutate(gene = str_match(gene, "g[0-9][0-9]?[0-9]?[0-9]?[0-9]?")) -> cazy_05_gene

#Drop the ".hmm" and paste all unique CAZY terms per gene
cazy_05_gene$domain_name <- sub(".hmm", "", cazy_05_gene$domain_name)
cazy_05_gene %>% 
  group_by(gene) %>% 
  distinct(domain_name, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  mutate(HMM = paste0(domain_name, collapse = ",")) %>% 
  select(gene, HMM) %>% 
  distinct(gene, .keep_all = TRUE) -> cazy_HMM

#make output table fs
cazy_05_gene %>% 
  group_by(gene) %>% 
  mutate(min = min(sequence_evalue)) %>% 
  filter(sequence_evalue == min) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  select(gene, domain_name) -> cazy_HMM_min

fs_cazy <- rownames_to_column(fs, var = "gene")
right_join(cazy_HMM_min, fs_cazy, by = "gene") -> fs_cazy
fs_cazy$anno <- paste(fs_cazy$gene, ":", fs_cazy$domain_name)

png(filename = "cellulases.png")
ggmaplot(fs_cazy, main = expression("Group 1" %->% "Group 2"),
         fdr = 0.05, fc = 2, size = 1, alpha = 0.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(fs_cazy$anno),
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         top = 0,
         label.select = c("g1056 : PL9_2", "g5926 : PL9_2", "g13764 : GH5_5", "g5517 : GH5_10", "g1703 : GH5_5", "g9411 : GH5_5"))
dev.off()
ggsave(filename = "cellulases.pdf", device = cairo_pdf())


#custom cellulases
sig <- fs_cazy$anno %in% c("g1056 : PL9_2", "g5926 : PL9_2", "g13764 : GH5_5", "g5517 : GH5_10", "g1703 : GH5_5", "g9411 : GH5_5")
genenames <- fs_cazy$gene
data <- data.frame(name = fs_cazy$anno, mean = fs_cazy$baseMean, lfc = fs_cazy$log2FoldChange,
                   padj = fs_cazy$padj, sig = sig)
data$mean <- log2(data$mean +1)
fc = 2
labs_data <- subset(data, sig)


p <- ggplot(data, aes(x = mean, y = lfc)) +
  geom_point(size = 1, alpha = 0.5, color = "darkgrey") +
  geom_point(data = filter(data, lfc > 1 & padj < 0.05), size = 1, color = "#87899E") +
  geom_point(data = filter(data, lfc > 2 & padj < 0.05), size = 1, color = "#5C5E6E") +
  scale_x_continuous(breaks=seq(0, max(data$mean), 2)) +
  labs(x = "Log2 mean expression", y = "Log2 fold change") +
  geom_hline(yintercept = c(0,  log2(fc)), linetype = c(1, 2),
             color = c("black", "black")) +
  geom_hline(yintercept = c(0,  log2(fcc)), linetype = c(1, 2),
             color = c("black",  "black")) +
  ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = as.character(name)),
                            box.padding = unit(0.25, "lines"),
                            point.padding = unit(0.5, "lines"),
                            force = 200, fontface = "bold",
                            size = 8, color = "black", label.padding = 0.5) +
  geom_point(data = labs_data, color = "red", size = 2) +
  theme_classic(base_size = 18)

#ggsave(filename = "cellulases.pdf", p)
ggsave(filename = "cellulases.png", p, device = "png", width = 9, height = 8, units = "in")

#custom cellulases - color
sig <- fs_cazy$anno %in% c("g1056 : PL9_2", "g5926 : PL9_2", "g13764 : GH5_5", "g5517 : GH5_10", "g1703 : GH5_5", "g9411 : GH5_5")
genenames <- fs_cazy$gene
data <- data.frame(name = fs_cazy$anno, mean = fs_cazy$baseMean, lfc = fs_cazy$log2FoldChange,
                   padj = fs_cazy$padj, sig = sig)
data$mean <- log2(data$mean +1)
fc = 2
fcc = 4
labs_data <- subset(data, sig)


p <- ggplot(data, aes(x = mean, y = lfc)) +
  geom_point(size = 1, alpha = 0.5, color = "darkgrey") +
  geom_point(data = filter(data, lfc > 1 & padj < 0.05), size = 1, color = "#87899E") +
  geom_point(data = filter(data, lfc > 2 & padj < 0.05), size = 1, color = "#5C5E6E") +
  scale_x_continuous(breaks=seq(0, max(data$mean), 2)) +
  labs(x = "Log2 mean expression", y = "Log2 fold change") +
  geom_hline(yintercept = c(0,  log2(fc)), linetype = c(1, 2),
             color = c("black", "black")) +
  geom_hline(yintercept = c(0,  log2(fcc)), linetype = c(1,  2),
             color = c("black",  "black")) +
  geom_point(data = labs_data, color = "red", size = 2) +
  theme_classic(base_size = 18)

ggsave(filename = "cellulases_nolabel.png", p, device = "png", width = 9, height = 8, units = "in")


#custom chitin
sig <- fs_cazy$anno %in% c("g4820 : CBM18", "g1183 : CBM50", "g5125 : AA11", "g9951 : CBM50", "g3373 : CBM18", "g1388 : GH18", "g2605 : CBM18")
genenames <- fs_cazy$gene
data <- data.frame(name = fs_cazy$anno, mean = fs_cazy$baseMean, lfc = fs_cazy$log2FoldChange,
                   padj = fs_cazy$padj, sig = sig)
data$mean <- log2(data$mean +1)
fc = 2
labs_data <- subset(data, sig)

p <- ggplot(data, aes(x = mean, y = lfc)) +
  geom_point(size = 1, alpha = 0.5, color = "darkgrey") +
  geom_point(data = filter(data, lfc > 1 & padj < 0.05), size = 1, color = "#87899E") +
  geom_point(data = filter(data, lfc > 2 & padj < 0.05), size = 1, color = "#5C5E6E") +
  geom_point(data = labs_data, color = "red", size = 2) +
  scale_x_continuous(breaks=seq(0, max(data$mean), 2)) +
  labs(x = "Log2 mean expression", y = "Log2 fold change") +
  geom_hline(yintercept = c(0,  log2(fc)), linetype = c(1,  2),
             color = c("black", "black")) +
  geom_hline(yintercept = c(0,  log2(fcc)), linetype = c(1,  2),
             color = c("black", "black")) +
  ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = as.character(name)),
                            box.padding = unit(0.25, "lines"),
                            point.padding = unit(0.5, "lines"),
                            force = 100, fontface = "bold",
                            size = 8, color = "black", label.padding = 0.5) +
  theme_classic(base_size = 18)

#(filename = "chitinases.pdf", p)
ggsave(filename = "chitinases.png", p, device = "png", width = 9, height = 8, units = "in")



#Test GH5_5
#res_feeding["g1703",]
#plotCounts(dds, "g1703")

#Extracting transformed count values
vsd <- vst(dds, blind = FALSE)

setwd("X:/210129_DGE_logFC1_DESeq2/PCA")
#PCA plot 
pdf(file="PCA.pdf")
plotPCA(vsd, intgroup=c("condition")) +
  ggtitle("PCA plot")
dev.off()

#PCA custom
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
levels(pcaData$condition) <- c("seeking", "feeding", "digesting")
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file="PCA.pdf", width = 5, height = 7)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=c("#3333CC", "#D60093", "#1F4E79")) +
  theme_minimal() +
  theme(panel.border = element_rect(fill = "NA"))
dev.off()

#Glimma plots without annotation
#fs
status <- as.numeric(res_fs$padj < .05 & abs(res_fs$log2FoldChange) > 1)
anno <- data.frame(GeneID=rownames(res_fs))
glMDPlot(res_fs, status = status, counts=(log(counts(dds,normalized=TRUE)+0.5)), groups=dds$condition, transform=FALSE, samples=colnames(dds),
         anno = anno, path=("~/transfer/Jenny/210129_DGE_logFC1_DESeq2"), folder="glimma-MD_fs", launch=FALSE)

#Glimma plots
#df
status <- as.numeric(res_df$padj < .05 & abs(res_df$log2FoldChange) > 1)
anno <- data.frame(GeneID=rownames(res_df))
glMDPlot(res_df, status = status, counts=(log(counts(dds,normalized=TRUE)+0.5)), groups=dds$condition, transform=FALSE, samples=colnames(dds),
         anno = anno, path=("~/transfer/Jenny/210129_DGE_logFC1_DESeq2"), folder="glimma-MD_df", launch=FALSE)

#Glimma plots
#ds
status <- as.numeric(res_ds$padj < .05 & abs(res_ds$log2FoldChange) > 1)
anno <- data.frame(GeneID=rownames(res_ds))
glMDPlot(res_ds, status = status, counts=(log(counts(dds,normalized=TRUE)+0.5)), groups=dds$condition, transform=FALSE, samples=colnames(dds),
         anno = anno, path=("~/transfer/Jenny/210129_DGE_logFC1_DESeq2"), folder="glimma-MD_ds", launch=FALSE)

#Define DGE gene set log FC 1
fs <- data.frame(res_fs)
df <- data.frame(res_df)
ds <- data.frame(res_ds)
fs %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) -> fs_1
df %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) -> df_1
ds %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) -> ds_1
dge_set <- c(rownames(fs_1), rownames(df_1), rownames(ds_1))
dge_set <- unique(dge_set)
#Generate z-score matrix for heatmaps
vsd_matrix <- assay(vsd)
vsd_matrix_dge <- vsd_matrix[dge_set,]
heat <- t(scale(t(vsd_matrix_dge)))
#heatmap(heat)

#Define DGE gene set log FC 2
fs %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 2) -> fs_2
df %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 2) -> df_2
ds %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 2) -> ds_2
dge_set_2 <- c(rownames(fs_2), rownames(df_2), rownames(ds_2))
dge_set_2 <- unique(dge_set_2)
#Generate z-score matrix for heatmaps
vsd_matrix_dge_2 <- vsd_matrix[dge_set_2,]
heat_2 <- t(scale(t(vsd_matrix_dge_2)))
#heatmap(heat_2)
#heatmap(vsd_matrix_dge_2)

#Clustering of z-scores with euclidian distance yields same clusters as clustering vst with Pearson
#correlation but at different height 
#par(mfrow = c(2,2))
#plot(hclust(dist(heat_2, method = "euclidean")))
#plot(hclust(dist(vsd_matrix_dge_2, method = "euclidean")))
#plot(hclust(as.dist(1 - cor(t(vsd_matrix_dge_2)))))


# Correlate like here: http://www.opiniomics.org/you-probably-dont-understand-heatmaps/
#Clusters rows by Pearson correlation as distance method.
hc <- hclust(as.dist(1 - cor(t(as.matrix(vsd_matrix_dge_2)))))
my_gene_partition_assignments2 <- cutree(hc, h=60/100*max(hc$height))

#make factor-labeling.txt for GO-seq
factor_labeling <- data.frame(my_gene_partition_assignments, dge_set_2)
factor_labeling$my_gene_partition_assignments <- paste0("c", factor_labeling$my_gene_partition_assignments)
setwd("X:/210129_DGE_logFC1_DESeq2/GO")
write.table(factor_labeling, file = "factor_labeling.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


#Generate subclusters for logFC 2
#? Generate them from normalised counts -> z-scores from the FC2 cutoff heat_2

#my_gene_dist <- dist(heat_2, method = "euclidean")
#my_hc_genes <- hclust(my_gene_dist, method = "complete")
#my_gene_partition_assignments <- cutree(my_hc_genes, h=80/100*max(my_hc_genes$height))


#Visualise clusters
'hc <- hclust(my_gene_dist, method = "complete")'
#hc_1 <- cutree(hc, h = 1)
hc_7.5 <- cutree(hc, h = 60/100*max(hc$height))
#hc_10 <- cutree(hc, h = 10)
the_bars <- cbind(hc_7.5)

setwd("X:/Talks/DGP2021/Figures")
png(filename = "clustering.png", width = 10, height = 5, unit = "in", res = 300)
plot(as.dendrogram(hc), leaflab = "none")
colored_bars(the_bars, as.dendrogram(hc), sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("h=60% height"),cex.rowLabels=0.7)
abline(h=60/100*max(hc$height), lty = 2, col="grey")
dev.off()

#make list
clusterlist <- list()

for (i in c(1:max(my_gene_partition_assignments))) {
  cluster <- heat_2[(my_gene_partition_assignments == i),]
  clusterlist[[i]] <- cluster
}

#Make line plots
setwd("X:/210129_DGE_logFC1_DESeq2/test")
pdf(file = "clusterplots.pdf")
par(mfrow=c(2, 2))
par(cex=0.6)
par(mar=c(7,4,4,2))

for (i in 1:max(my_gene_partition_assignments)) {
  data = clusterlist[[i]]
  ymin = min(data); ymax = max(data);
  plot_label = paste("cluster", i, ', ', length(data[,1]), "genes", sep='')
  plot(as.numeric(data[1,]), type='l', main = plot_label, ylim=c(ymin,ymax), col='lightgray', xaxt='n', xlab='', ylab='z-score')
  axis(side=1, at=1:length(data[1,]), labels=colnames(data), las=2)
  for(r in 2:length(data[,1])) {
    points(as.numeric(data[r,]), type='l', col='lightgray')
  }
  points(as.numeric(colMeans(data)), type='l', col='darkblue')
}
dev.off()

#Make boxplots

p <- list()
for (i in 1:max(my_gene_partition_assignments)) {
  test <- rownames_to_column(as.data.frame(clusterlist[[i]]), var = "gene")
  test_longer <- pivot_longer(test, cols = !gene, names_to = "sample")
  test_longer %>% 
    filter(sample == "V1S4" | sample == "V1S5" | sample == "V1S6") %>% 
    mutate(cat = "Seeking") -> test_longer1
  test_longer %>% 
   filter(sample == "V1S7" | sample == "V1S8" | sample == "V1S9") %>% 
    mutate(cat = "Feeding") -> test_longer2
  test_longer %>% 
    filter(sample == "V1S1" | sample == "V1S2" | sample == "V1S3") %>% 
    mutate(cat = "Digesting") -> test_longer3
  comb <- rbind(test_longer1, test_longer2, test_longer3)
  comb$cat <-factor(comb$cat, levels = c("Seeking", "Feeding", "Digesting"))
  g <- ggplot(comb, aes(x = cat, y = value, fill = cat)) +
        geom_boxplot(outlier.colour = "darkgray",
                     outlier.size=0.1,
                     outlier.shape = 1,
                     outlier.alpha = 0.5,
                     outlier.stroke = 0.1,
                     size=0.1,
                     colour="black",
                     notch=TRUE) +
        ggtitle(paste("cluster", i, "|", "genes:", nrow(clusterlist[[i]]))) +
        ylab("z-score") +
        xlab(NULL) +
        theme(legend.position="none") +
        scale_fill_brewer(palette="Set2")
  p[[i]] <- ggplotGrob(g)
}
pdf(file = "boxplots.pdf", width = 15, height = 9)
grid.arrange(grobs = p, ncol=3)
dev.off()


#Make boxplots - DGEPatterns look

p <- list()
splan <- 3 - 1L

for (i in 1:max(my_gene_partition_assignments)) {
  test <- rownames_to_column(as.data.frame(clusterlist[[i]]), var = "gene")
  test_longer <- pivot_longer(test, cols = !gene, names_to = "sample")
  test_longer %>% 
    filter(sample == "V1S4" | sample == "V1S5" | sample == "V1S6") %>% 
    mutate(cat = "Seeking") -> test_longer1
  test_longer %>% 
    filter(sample == "V1S7" | sample == "V1S8" | sample == "V1S9") %>% 
    mutate(cat = "Attacking") -> test_longer2
  test_longer %>% 
    filter(sample == "V1S1" | sample == "V1S2" | sample == "V1S3") %>% 
    mutate(cat = "Digesting") -> test_longer3
  comb <- rbind(test_longer1, test_longer2, test_longer3)
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
pdf(file = "blue_boxplots.pdf", width = 15, height = 3)
grid.arrange(grobs = p, ncol=5)
dev.off()

pdf(file = "blue_boxplots_rowwise.pdf", width = 3, height = 15)
grid.arrange(grobs = p, nrow=5)
dev.off()

#violin

for (i in 1:max(my_gene_partition_assignments)) {
  test <- rownames_to_column(as.data.frame(clusterlist[[i]]), var = "gene")
  test_longer <- pivot_longer(test, cols = !gene, names_to = "sample")
  test_longer %>% 
    filter(sample == "V1S4" | sample == "V1S5" | sample == "V1S6") %>% 
    mutate(cat = "Seeking") -> test_longer1
  test_longer %>% 
    filter(sample == "V1S7" | sample == "V1S8" | sample == "V1S9") %>% 
    mutate(cat = "Feeding") -> test_longer2
  test_longer %>% 
    filter(sample == "V1S1" | sample == "V1S2" | sample == "V1S3") %>% 
    mutate(cat = "Digesting") -> test_longer3
  comb <- rbind(test_longer1, test_longer2, test_longer3)
  g <- ggplot(comb, aes(x = cat, y = value, fill = "#F8766D", color = "#F8766D")) +
    geom_violin(outlier.size = 0,
                 outlier.shape = NA,
                 alpha = 0.5) +
    geom_point(alpha = 0.2, size = .5,
               position = position_jitterdodge(dodge.width = 0.9)) +
    stat_smooth(aes(x = cat, y = value,
                    group = "#F8766D", color = "#F8766D"),
                se = TRUE,
                method = "lm", formula = y~poly(x, splan)) +
    geom_line(alpha = 0.2) +
    ggtitle(paste("cluster", i, "|", "genes:", nrow(clusterlist[[i]]))) +
    ylab("z-score") +
    xlab(NULL) +
    theme(legend.position="none")
  p[[i]] <- ggplotGrob(g)
}
pdf(file = "violinplots.pdf", width = 15, height = 9)
grid.arrange(grobs = p, ncol=3)
dev.off()



#Continue on local PC with ComplexHeatmap

#Parse IPS annotation for later
setwd("C:/Users/Guest/Desktop")
ips <- read.delim("IPS_isoforms.tsv", header = TRUE, col.names = c("qseqid" , "md5" , "length" , "analysis", "accession", "description", "start", "stop", "evalue", "status", "date","IPS_accession", "IPS_description", "GO", "PA"))
ips$gene <- sub("_i[0-9][0-9]?[0-9]?.p.", "", ips$qseqid)
ips %>% 
  mutate(gene = str_match(gene, "g[0-9][0-9]?[0-9]?[0-9]?[0-9]?")) %>%
  select(gene, description) %>% 
  filter(description != "-") %>% 
  group_by(gene) %>% 
  distinct(description, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  mutate(ips_description = paste0(description, collapse = ",")) %>% 
  group_by(gene) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  select(gene, description) -> ips_pergene

nrow(ips_pergene)
remove(ips)
#Max character count 350
#ips_pergene$ips_description <- strtrim(ips_pergene$ips_description, 350)

#Read in and parse eggnog annotation
#HMM
setwd("X:/210127_eggnog_HMMER")
eggnog_HMM <- read.delim("eggnog_HMMERout.tab_EggNOGannot_evalLessThan1e-10_noheader.txt", header = TRUE, sep = "\t")
#Create gene column 
eggnog_HMM$gene <- sub("_i[0-9][0-9]?[0-9]?.p.", "", eggnog_HMM$SequenceName)
eggnog_HMM %>% 
  mutate(gene = str_match(gene, "g[0-9][0-9]?[0-9]?[0-9]?[0-9]?")) -> eggnog_HMM_gene 

#Diamond
setwd("X:/210205_another_eggnog")
eggnog_diamond <- read.delim("eggnog_5.0_diamond.txt", header = TRUE, sep = "\t")

#Create gene column
eggnog_diamond$gene <- sub("_i[0-9][0-9]?[0-9]?.p.", "", eggnog_diamond$X.query_name)
eggnog_diamond %>%
  mutate(gene = str_match(gene, "g[0-9][0-9]?[0-9]?[0-9]?[0-9]?")) -> eggnog_diamond_gene 

eggnog_HMM_gene %>% 
  select(gene, COG, description) %>% 
  group_by(gene) %>% 
  distinct(COG, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  mutate(COG2 = paste0(COG, collapse = "")) %>% 
  group_by(gene) %>% 
  mutate(description2 = paste0(description, collapse = "")) %>% 
  select(gene, COG2, description2) %>% 
  group_by(gene) %>% 
  distinct(COG2, .keep_all = TRUE) -> eggnog_HMM_select

eggnog_diamond_gene %>% 
  select(gene, seed_eggNOG_ortholog, Preferred_name, best_og_cat, best_og_desc) %>% 
  filter(best_og_cat != "-") %>% 
  group_by(gene) %>% 
  distinct(best_og_cat, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  mutate(COG = paste0(best_og_cat, collapse = "")) %>% 
  group_by(gene) %>% 
  mutate(description = paste0(best_og_desc, collapse = "")) %>% 
  select(gene, seed_eggNOG_ortholog, Preferred_name, COG, description) %>% 
  group_by(gene) %>% 
  distinct(COG, .keep_all = TRUE) -> eggnog_diamond_select


sum(eggnog_HMM_select$gene %in% eggnog_diamond_select$gene)
#7465 / 7567
sum(eggnog_diamond_select$gene %in% eggnog_HMM_select$gene)
#7465 / 11415

#is in Diamond?
#notkeep <- eggnog_HMM_select$gene %in% eggnog_diamond_select$gene
#eggnog_HMM_select <- eggnog_HMM_select[!notkeep,]

#eggnog_merged <- rbind(eggnog_HMM_select, eggnog_diamond_select)
#eggnog_merged %>% 
  filter(COG2 != "") -> eggnog_merged 
eggnog_merged <- eggnog_diamond_select
nrow(eggnog_merged[eggnog_merged$Preferred_name != "-",])/nrow(eggnog_merged)
  
heat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  left_join(eggnog_merged) %>% 
  filter(COG != "") %>% 
  mutate(anno = paste(gene, " | ", Preferred_name, " | ", description)) -> heat_anno

mylist <- list()
for (i in LETTERS[-24]){
  a <- filter(heat_anno, str_detect(COG, i))
  mylist[[length(mylist) + 1]] <- a
}
names(mylist) <- LETTERS[-24]

setwd("X:/210129_DGE_logFC1_DESeq2")

eggnog_terms = c("A: RNA processing and modification", "B: Chromatin structure and dynamics", "C: Energy production and conversion", "D: Cell cycle control, cell division, chromosome partitioning", "E: Amino acid transport and metabolism", "F: Nucleotide transport and metabolism", "G: Carbohydrate transport and metabolism", "H: Coenzyme transport and metabolism", "I: Lipid transport and metabolism", "J: Translation, ribosomal structure and biogenesis", "K: Transcription", "L: Replication, recombination and repair", "M: Cell wall/membrane/envelope biogenesis", "N: Cell motility", "O: Posttranslational modification, protein turnover, chaperones", "P: Inorganic ion transport and metabolism", "Q: Secondary metabolites biosynthesis, transport and catabolism", "R: General function prediction only", "S: Function unknown", "T: Signal transduction mechanisms", "U: Intracellular trafficking, secretion, and vesicular transport", "V: Defense mechanisms", "W: Extracellular structures", "Y: Nuclear structure", "Z: Cytoskeleton")



#Create one heatmap per COG functional cluster
for (i in c(1:25)){
  mylist[[i]] %>% 
    select(c(15,8,9,10,2,3,4,5,6,7)) %>% 
    column_to_rownames(var = "anno") -> b
  if (nrow(mylist[[i]] > 0)){
    plot_height <- nrow(mylist[[i]])/4
    if (plot_height < 4){
      plot_height <- 4
    } else {
      plot_height <- nrow(mylist[[i]])/4
    }
    pdf(file = gsub(" ", "", paste(names(mylist)[i], ".pdf")), height = plot_height, width = 24)
    draw(Heatmap(as.matrix(b), 
                 name = "z-score", 
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 column_title = c(paste(eggnog_terms[i], paste("Number of DGE: ", nrow(mylist[[i]])), sep="\n")), 
                 column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                 show_row_dend = FALSE,
                 column_split = rep(1:3, each = 3),
                 bottom_annotation = HeatmapAnnotation(foo = anno_block(labels = c("Digesting", "Starving", "Feeding"))),
                 col = colorRamp2(c(-2,-1,0,1,2), rev(brewer.pal(5, "RdYlBu")))), 
         padding = unit(c(2,2,2,400), "mm"), 
         heatmap_legend_side = "left")
    dev.off()
  } else {
    next
  }
}

#Make boxplot of z-scores per functional cluster
p_COG <- list()
for (i in c(1:25)) {
  mylist[[i]] %>% 
    select(c(8,9,10,2,3,4,5,6,7)) -> c
    c <- rownames_to_column(as.data.frame(c), var = "gene")
    c_longer <- pivot_longer(c, cols = !gene, names_to = "sample")
    c_longer %>% 
      filter(sample == "V1S4" | sample == "V1S5" | sample == "V1S6") %>% 
      mutate(cat = "Seeking") -> c_longer1
    c_longer %>% 
      filter(sample == "V1S7" | sample == "V1S8" | sample == "V1S9") %>% 
      mutate(cat = "Feeding") -> c_longer2
    c_longer %>% 
      filter(sample == "V1S1" | sample == "V1S2" | sample == "V1S3") %>% 
      mutate(cat = "Digesting") -> c_longer3
    ccomb <- rbind(c_longer1, c_longer2, c_longer3)
    g <- ggplot(ccomb, aes(x = cat, y = value, fill = cat)) +
      geom_boxplot(outlier.colour = "darkgray",
                   outlier.size=0.1,
                   outlier.shape = 1,
                   outlier.alpha = 0.5,
                   outlier.stroke = 0.1,
                   size=0.1,
                   colour="black") +#,
      #notch=TRUE) +
      ggtitle(paste(eggnog_terms[i], "|", "genes:", nrow(mylist[[i]]))) +
      ylab("z-score") +
      xlab(NULL) +
      theme(legend.position="none") +
      scale_fill_brewer(palette="Set2")
    p_COG[[i]] <- ggplotGrob(g)
}

pdf(file = "COG_boxplots.pdf", width = 20, height = 20)
grid.arrange(grobs = p_COG[lengths(p_COG) != 0])
dev.off()

#make heatmap with COG mean z-scores

colmeans <- list()
for (i in c(1:25)) {
  colmeans[[i]] <- colMeans(select(mylist[[i]], c(2:10)))
}
colmeans_df <- do.call(rbind, colmeans)
rownames(colmeans_df) <- eggnog_terms

pdf("COG_mean_z.pdf", width = 10, height = 8)
draw(Heatmap(as.matrix(colmeans_df[-c(18),]), 
             name = "mean z-score", 
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_column_names = FALSE,
             column_split = rep(1:3, each = 3),
             bottom_annotation = HeatmapAnnotation(foo = anno_block(labels = c("Seeking", "Feeding", "Digesting"))),
             col = colorRamp2(c(-1,-0.5,0,0.5,1), rev(brewer.pal(5, "RdYlBu")))),
     padding = unit(c(2,2,2,40), "mm"), 
     heatmap_legend_side = "left")
dev.off()

#Import CAZY annotation
setwd("X:/210129_DGE_logFC1_DESeq2/Annotation")
cazy <- read_tblout("CAZY_noEval.out")
cazy %>% 
  filter(best_domain_evalue < 1E-4) %>% 
  dplyr::select(domain_name, query_name, sequence_evalue) -> cazy_05

#Add gene columns
cazy_05$gene <- sub("_i[0-9][0-9]?[0-9]?.p.", "", cazy_05$query_name)
cazy_05 %>% 
  mutate(gene = str_match(gene, "g[0-9][0-9]?[0-9]?[0-9]?[0-9]?")) -> cazy_05_gene

#Drop the ".hmm" and paste all unique CAZY terms per gene
cazy_05_gene$domain_name <- sub(".hmm", "", cazy_05_gene$domain_name)
cazy_05_gene %>% 
  group_by(gene) %>% 
  distinct(domain_name, .keep_all = TRUE) %>% 
  group_by(gene) %>% 
  mutate(HMM = paste0(domain_name, collapse = ",")) %>% 
  select(gene, HMM) %>% 
  distinct(gene, .keep_all = TRUE) -> cazy_HMM

#make output table fs
cazy_05_gene %>% 
  group_by(gene) %>% 
  mutate(min = min(sequence_evalue)) %>% 
  filter(sequence_evalue == min) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  select(gene, domain_name) -> cazy_HMM_min

fs_cazy <- rownames_to_column(fs, var = "gene")
right_join(cazy_HMM_min, fs_cazy, by = "gene") -> fs_cazy
fs_cazy$anno <- paste(fs_cazy$gene, ":", fs_cazy$domain_name)
write.table(fs_cazy, file = "fs_cazy", quote = FALSE, sep = "\t")

#Generate input for HMM
heat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  inner_join(cazy_HMM) %>% 
  mutate(anno = paste(gene, " | ", HMM)) -> heat_cazy

#Make list for heatmaps
cazy_groups <- c("GH", "GT", "PL", "CE", "AA", "CBM")
mylist_cazy <- list()
for (i in cazy_groups){
  a <- filter(heat_cazy, str_detect(HMM, i))
  mylist_cazy[[length(mylist_cazy) + 1]] <- a
}
names(mylist_cazy) <- cazy_groups

#Create one heatmap per CAZY group

cazy_title <- c("Glycoside Hydrolases (GHs)", " Glycosyl Transferases (GTs)", "Polysaccharide Lyases (PLs)", "Carbohydrate Esterases (CEs)", "Auxiliary Activities (AAs)", "Carbohydrate-Binding Modules (CBMs)")

for (i in c(1:6)){
  mylist_cazy[[i]] %>% 
    select(c(12,8,9,10,2,3,4,5,6,7)) %>% 
    column_to_rownames(var = "anno") -> b
  if (nrow(mylist_cazy[[i]] > 0)){
    plot_height <- nrow(mylist_cazy[[i]])/4
    if (plot_height < 4){
      plot_height <- 2
    } else {
      plot_height <- nrow(mylist_cazy[[i]])/4
    }
    pdf(file = gsub(" ", "", paste(names(mylist_cazy)[i], ".pdf")), height = plot_height, width = 24)
    draw(Heatmap(as.matrix(b), 
                 name = "z-score", 
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 column_title = c(paste(cazy_title[i], paste("Number of DGE: ", nrow(mylist_cazy[[i]])), sep="\n")), 
                 column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                 show_row_dend = FALSE,
                 column_split = rep(1:3, each = 3),
                 bottom_annotation = HeatmapAnnotation(foo = anno_block(labels = c("Digesting", "Starving", "Feeding"))),
                 col = colorRamp2(c(-2,-1,0,1,2), rev(brewer.pal(5, "RdYlBu")))), 
         padding = unit(c(2,2,2,400), "mm"), 
         heatmap_legend_side = "left")
    dev.off()
  } else {
    next
  }
}

#Glimma plots without annotation
setwd("X:/210129_DGE_logFC1_DESeq2/Glimma_anno")
#fs
status <- as.numeric(res_fs$padj < .05 & abs(res_fs$log2FoldChange) > 1)
anno <- data.frame(gene = rownames(res_fs))
anno %>% 
  left_join(eggnog_merged, by = "gene") %>% 
  left_join(ips_pergene, by = "gene") %>% 
  left_join(cazy_HMM, by = "gene") -> anno
colnames(anno) <- c("GeneID", "COG", "eggnog_description", "ips_description", "CAZy")
glMDPlot(res_fs, status = status, counts=(log(counts(dds,normalized=TRUE)+0.5)), groups=dds$condition, transform=FALSE, samples=colnames(dds),
         anno = anno, path=("X:/210129_DGE_logFC1_DESeq2/Glimma_anno"), folder="glimma-MD_fs", launch=FALSE)


#df
status <- as.numeric(res_df$padj < .05 & abs(res_df$log2FoldChange) > 1)
anno <- data.frame(gene = rownames(res_df))
anno %>% 
  left_join(eggnog_merged, by = "gene") %>% 
  left_join(ips_pergene, by = "gene") %>% 
  left_join(cazy_HMM, by = "gene") -> anno
colnames(anno) <- c("GeneID", "COG", "eggnog_description", "ips_description", "CAZy")
glMDPlot(res_df, status = status, counts=(log(counts(dds,normalized=TRUE)+0.5)), groups=dds$condition, transform=FALSE, samples=colnames(dds),
         anno = anno, path=("X:/210129_DGE_logFC1_DESeq2/Glimma_anno"), folder="glimma-MD_df", launch=FALSE)

#ds
status <- as.numeric(res_ds$padj < .05 & abs(res_ds$log2FoldChange) > 1)
anno <- data.frame(gene = rownames(res_ds))
anno %>% 
  left_join(eggnog_merged, by = "gene") %>% 
  left_join(ips_pergene, by = "gene") %>% 
  left_join(cazy_HMM, by = "gene") -> anno
colnames(anno) <- c("GeneID", "COG", "eggnog_description", "ips_description", "CAZy")
glMDPlot(res_ds, status = status, counts=(log(counts(dds,normalized=TRUE)+0.5)), groups=dds$condition, transform=FALSE, samples=colnames(dds),
         anno = anno, path=("X:/210129_DGE_logFC1_DESeq2/Glimma_anno"), folder="glimma-MD_ds", launch=FALSE)


#Try out with normalised counts
#Make boxplot of vst per functional cluster
setwd("X:/210129_DGE_logFC1_DESeq2/test/vst")
#make input

vsd_matrix_dge %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  left_join(eggnog_merged) %>% 
  filter(COG != "") %>% 
  mutate(anno = paste(gene, " | ", Preferred_name, " | ", description)) -> heat_vst

mylist_vst <- list()
for (i in LETTERS[-24]){
  a <- filter(heat_vst, str_detect(COG, i))
  mylist_vst[[length(mylist_vst) + 1]] <- a
}
names(mylist_vst) <- LETTERS[-24]

colnames(mylist_vst[[1]])

pp_COG <- list()
for (i in c(1:25)) {
  mylist_vst[[i]] %>% 
    select(c(8,9,10,2,3,4,5,6,7)) -> c
  c <- rownames_to_column(as.data.frame(c), var = "gene")
  c_longer <- pivot_longer(c, cols = !gene, names_to = "sample")
  c_longer %>% 
    filter(sample == "V1S4" | sample == "V1S5" | sample == "V1S6") %>% 
    mutate(cat = "Seeking") -> c_longer1
  c_longer %>% 
    filter(sample == "V1S7" | sample == "V1S8" | sample == "V1S9") %>% 
    mutate(cat = "Feeding") -> c_longer2
  c_longer %>% 
    filter(sample == "V1S1" | sample == "V1S2" | sample == "V1S3") %>% 
    mutate(cat = "Digesting") -> c_longer3
  ccomb <- rbind(c_longer1, c_longer2, c_longer3)
  g <- ggplot(ccomb, aes(x = cat, y = value, fill = cat)) +
    geom_boxplot(outlier.colour = "darkgray",
                 outlier.size=0.1,
                 outlier.shape = 1,
                 outlier.alpha = 0.5,
                 outlier.stroke = 0.1,
                 size=0.1,
                 colour="black") +#,
    #notch=TRUE) +
    ggtitle(paste(eggnog_terms[i], "|", "genes:", nrow(mylist_vst[[i]]))) +
    ylab("z-score") +
    xlab(NULL) +
    theme(legend.position="none") +
    scale_fill_brewer(palette="Set2")
  pp_COG[[i]] <- ggplotGrob(g)
}

pdf(file = "COG_boxplots.pdf", width = 20, height = 20)
grid.arrange(grobs = pp_COG[lengths(p_COG) != 0])
dev.off()



'#Single plot for testing purposes
mylist$Z %>% 
  select(c(13,8,9,10,2,3,4,5,6,7)) %>% 
  column_to_rownames(var = "anno") -> b
draw(Heatmap(as.matrix(b), 
             name = "z-score", 
             cluster_columns = FALSE,
             show_column_names = FALSE,
             column_title = NULL,
             show_row_dend = FALSE,
             column_split = rep(1:3, each = 3),
             bottom_annotation = HeatmapAnnotation(foo = anno_block(labels = c("Digesting", "Starving", "Feeding"))),
             col = colorRamp2(c(-2,-1,0,1,2), rev(brewer.pal(5, "RdYlBu")))), 
     padding = unit(c(2,2,2,400), "mm"), 
     heatmap_legend_side = "left")'


'#Works, two line breaks in column title
for (i in c(1:26)){
  mylist[[i]] %>% 
    select(c(13,8,9,10,2,3,4,5,6,7)) %>% 
    column_to_rownames(var = "anno") -> b
  if (nrow(mylist[[i]] > 0)){
    plot_height <- nrow(mylist[[i]])/4
    pdf(file = gsub(" ", "", paste(names(mylist)[i], ".pdf")), height = plot_height, width = 24)
    draw(Heatmap(as.matrix(b), 
                 name = "z-score", 
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 column_title = c(paste(eggnog_terms[i], "Number of DGE: ", nrow(mylist[[i]]), sep="\n")), 
                 column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                 show_row_dend = FALSE,
                 column_split = rep(1:3, each = 3),
                 bottom_annotation = HeatmapAnnotation(foo = anno_block(labels = c("Digesting", "Starving", "Feeding"))),
                 col = colorRamp2(c(-2,-1,0,1,2), rev(brewer.pal(5, "RdYlBu")))),
         padding = unit(c(2,2,2,400), "mm"), 
         heatmap_legend_side = "left")
    dev.off()
  } else {
    next
  }
}'

'for (i in c(1:26)){
  mylist[[i]] %>% 
    select(c(13,8,9,10,2,3,4,5,6,7)) %>% 
    column_to_rownames(var = "anno") -> b
  if (nrow(mylist[[i]] > 0)){
      pdf(file = paste(i, ".pdf"))
      heatmap(as.matrix(b))
      dev.off()
  } else {
    next
  }
}'

'for (i in c(1:26)){
  mylist[[i]] %>% 
    select(c(13,8,9,10,2,3,4,5,6,7)) %>% 
    column_to_rownames(var = "anno") -> b
  if (nrow(mylist[[i]] > 0)){
    pdf(file = gsub(" ", "", paste(names(mylist)[i], ".pdf")), height = 24, width = 24)
    heatmap(as.matrix(b), Colv = FALSE)
    dev.off()
  } else {
    next
  }
}'
