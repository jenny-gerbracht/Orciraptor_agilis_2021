library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(goseq)
library(GO.db)
library(qvalue)
library(clusterProfiler)
library(forcats)

source(file = "../../config.txt")

###################
## Read GO annotation from Blast2GO
###################

blast2go <- read.delim(file = "blast2go_annot.annot",
                       header = FALSE,
                       sep = "\t",
                       col.names = c("peptide", "GO", "text"))
blast2go %>% 
  dplyr::select(peptide, GO) %>% 
  filter(str_detect(GO, "GO:")) -> blast2go
blast2go$gene <- str_split(blast2go$peptide, "_", simplify = TRUE)[,1]
#blast2go %>% 
#  group_by(gene) %>% 
#  summarise(all_GO = paste0(GO, collapse = ",")) -> go_annotation

blast2go %>% 
  dplyr::select(gene, GO) -> blast2go

###################
## Run Trinity scripts fasta_seq_length.pl and TPM_weighted_gene_length.py to obtain gene lengths
## Run Trinity script run_GOseq.pl to obtain goseq script
###################

###################
## Modified _runGOseq.R script from un_GOseq.pl
###################

# capture list of genes for functional enrichment testing
factor_labeling = read.delim(file = paste0(mydir, "/Module_5/GO_analysis/factor_labeling.txt"),
                             col.names = c("gene", "cluster"),
                             sep = "\t",
                             quote = "")
DE_genes = rownames(factor_labeling)


# get gene lengths
gene_lengths = read.table(file = paste0(mydir, "/Module_5/GO_analysis/gene_lengths.txt"), 
                          header = TRUE, 
                          row.names=1, 
                          com = '')
gene_lengths = as.matrix(gene_lengths[,1, drop=FALSE])


# get background gene list
background = read.table(file = paste0(mydir, "/Module_5/background.txt"),
                        row.names = 1)


# Cluster with cluster profiler
"cluster_dotplots <- list()
cluster_results <- list()
cluster_objects <- list()
clusters <- unique(factor_labeling$cluster)
sample_set_gene_ids = rownames(background)
sample_set_gene_lengths <- gene_lengths[sample_set_gene_ids,]"

"for (i in clusters) {
  deg <- filter(factor_labeling, cluster == (paste(i)))$gene
  cat_genes_vec = as.integer(sample_set_gene_ids %in% deg)
  pwf = nullp(cat_genes_vec, bias.data = sample_set_gene_lengths)
  rownames(pwf) = sample_set_gene_ids
  result <- goseq(pwf, gene2cat = blast2go, use_genes_without_cat = FALSE)
  result <- result[c(1:200),]
  #result[result$over_represented_pvalue > 1 - 1e-10,]$over_represented_pvalue <- c(1 - 1e-10)
  result$padj <- p.adjust(result$over_represented_pvalue, method="BH")

  GOs <- result$category
  geneSets <- list()

  for (j in GOs) {
    filter(go_annotation, str_detect(all_GO, paste(j))) %>% 
      summarise(genes = paste0(gene, collapse = "," )) -> tmp
    geneSets[[j]] <- str_split(tmp$genes, ",", simplify = TRUE)[str_split(tmp$genes, ",", simplify = TRUE) %in% filter(factor_labeling, cluster == paste(i))$gene]
  }

  results_df <- as.data.frame(matrix(vector(), nrow(result), 9,))
  colnames(results_df) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue", "geneID","Count")
  results_df$ID <- result$category
  results_df$Description <- result$term
  results_df$GeneRatio <- paste0(result$numDEInCat, "/", nrow(filter(factor_labeling, cluster == paste(i))))
  #results_df$GeneRatio <- result$numDEInCat/nrow(filter(factor_labeling, cluster == paste(i)))
  results_df$BgRatio <- paste0(result$numDEInCat, "/", nrow(filter(factor_labeling, cluster == paste(i))))
  results_df$pvalue <- result$over_represented_pvalue
  results_df$p.adjust <- p.adjust(result$over_represented_pvalue, method="BH")
  #results_df$qvalue <- qvalue(result$over_represented_pvalue)$qvalues 
  results_df$geneID <- geneSets[results_df$ID]
  results_df$Count <- result$numDEInCat
  
  my_object <- new("enrichResult",
                  readable = FALSE,
                  result = results_df,
                  pvalueCutoff = 1,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1,
                  organism = "UNKNOWN",
                  ontology = "UNKNOWN",
                  gene = filter(factor_labeling, cluster == paste(i))$gene ,
                  keytype = "UNKNOWN",
                  universe = background.gene_ids,
                  gene2Symbol = character(0),
                  geneSets = geneSets)
  
  cluster_objects[[i]] <- my_object
  cluster_results[[i]] <- results_df
  cluster_dotplots[[i]] <- enrichplot::dotplot(my_object, showCategory=30)
}
"
# Cluster with ggplot2
cluster_dotplots <- list()
cluster_results <- list()
clusters <- unique(factor_labeling$cluster)
sample_set_gene_ids = rownames(background)
sample_set_gene_lengths <- gene_lengths[sample_set_gene_ids,]

for (i in clusters) {
  deg <- filter(factor_labeling, cluster == i)$gene
  cat_genes_vec = as.integer(sample_set_gene_ids %in% deg)
  pwf = nullp(cat_genes_vec, bias.data = sample_set_gene_lengths)
  rownames(pwf) = sample_set_gene_ids
  result <- goseq(pwf, gene2cat = blast2go, use_genes_without_cat = FALSE)
  result <- result[c(1:200),]
  #result[result$over_represented_pvalue > 1 - 1e-10,]$over_represented_pvalue <- c(1 - 1e-10)
  result$padj <- p.adjust(result$over_represented_pvalue, method="BH")
  
  GOs <- result$category
  geneSets <- list()
  
  for (j in GOs) {
    filter(go_annotation, str_detect(all_GO, paste(j))) %>% 
      summarise(genes = paste0(gene, collapse = "," )) -> tmp
    geneSets[[j]] <- str_split(tmp$genes, ",", simplify = TRUE)[str_split(tmp$genes, ",", simplify = TRUE) %in% filter(factor_labeling, cluster == paste(i))$gene]
  }
  
  results_df <- as.data.frame(matrix(vector(), nrow(result), 9,))
  colnames(results_df) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue", "geneID","Count")
  results_df$ID <- result$category
  results_df$Description <- result$term
  #results_df$GeneRatio <- paste0(result$numDEInCat, "/", nrow(filter(factor_labeling, cluster == paste(i))))
  results_df$GeneRatio <- result$numDEInCat/nrow(filter(factor_labeling, cluster == i))
  results_df$BgRatio <- paste0(result$numDEInCat, "/", nrow(filter(factor_labeling, cluster == i)))
  results_df$pvalue <- result$over_represented_pvalue
  results_df$p.adjust <- p.adjust(result$over_represented_pvalue, method="BH")
  #results_df$qvalue <- qvalue(result$over_represented_pvalue)$qvalues 
  results_df$geneID <- geneSets[results_df$ID]
  results_df$Count <- result$numDEInCat
  
  
  cluster_results[[i]] <- results_df
  g <- ggplot(results_df[c(1:30),], aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
    geom_point(aes(size = Count, color = p.adjust)) +
    theme_bw(base_size = 14) +
    labs(size = "Count") +
    scale_colour_gradient(limits = c(0,0.05), low="red") +
    scale_size_continuous(limits = c(1, 160), range = c(2,10), breaks = c(1, 10, 50, 100, 150)) +
    ylab(NULL) +
    guides(size = guide_legend(order = 1)) +
    ggtitle(paste0("Cluster ", i))
  cluster_dotplots[[i]] <- ggplotGrob(g)
  
}

#pdf(file = paste0(mydir, "/Figures/cluster_dotplots_log2FC1.pdf"), width = 15, height = 30)
#grid.arrange(grobs = cluster_dotplots, nrow = 5)
#dev.off()

pdf(file = paste0(mydir, "/Figures/cluster_dotplots_log2FC1.pdf"), width = 14, height = 30)
grid::grid.newpage()
grid::grid.draw(rbind(cluster_dotplots[[1]], 
                      cluster_dotplots[[2]],
                      cluster_dotplots[[3]],
                      cluster_dotplots[[4]],
                      cluster_dotplots[[5]]))
dev.off()

View(cluster_results[[3]])


## Solution with ggplot 
## from Tommy's code

ggplot(cluster_results[[5]][c(1:30),], aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  labs(size = "Count") +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")


## count the gene number
gene_count<- x@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)

## merge with the original dataframe
dot_df<- left_join(x@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

## plot
library(forcats) ## for reordering the factor
ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")



## trinity

factor_labeling = read.delim(file = paste0(mydir, "/Module_5/GO_analysis/factor_labeling.txt"),
                             col.names = c("type"),
                             sep = "\t",
                             quote = "")
DE_genes = rownames(factor_labeling)


# get gene lengths
gene_lengths = read.table(file = paste0(mydir, "/Module_5/GO_analysis/gene_lengths.txt"), 
                          header = TRUE, 
                          row.names=1, 
                          com = '')
gene_lengths = as.matrix(gene_lengths[,1, drop=FALSE])


# get background gene list
background = read.table(file = paste0(mydir, "/Module_5/background.txt"),
                        row.names = 1)



# parse GO assignments
#GO_info = read.table("go_annotations.txt", header=F, row.names=1,stringsAsFactors=F)
GO_info_listed = apply(go_annotation, 1, function(x) unlist(strsplit(x,',')))
names(GO_info_listed) = go_annotation$gene
get_GO_term_descr =  function(x) {
  d = 'none';
  go_info = GOTERM[[x]];
  if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
  return(d);
}


#organize go_id -> list of genes
GO_to_gene_list = list()
for (gene_id in intersect(names(GO_info_listed), sample_set_gene_ids)) {
  go_list = GO_info_listed[[gene_id]]
  for (go_id in go_list) {
    GO_to_gene_list[[go_id]] = c(GO_to_gene_list[[go_id]], gene_id)
  }
}


# GO-Seq protocol: build pwf based on ALL DE features
missing_gene_lengths = sample_set_gene_ids[! sample_set_gene_ids %in% rownames(gene_lengths)]
if (length(missing_gene_lengths) > 0) {
  stop("Error, missing gene lengths for features: ", paste(missing_gene_lengths, collapse=', '))
}
sample_set_gene_lengths = gene_lengths[sample_set_gene_ids,]
GO_info_listed = GO_info_listed[ names(GO_info_listed) %in% sample_set_gene_ids ]
cat_genes_vec = as.integer(sample_set_gene_ids %in% rownames(factor_labeling))
pwf=nullp(cat_genes_vec, bias.data=sample_set_gene_lengths)
rownames(pwf) = sample_set_gene_ids


# perform functional enrichment testing for each category.

go_list <- list()

for (feature_cat in factor_list) {
  message('Processing category: ', feature_cat)
  gene_ids_in_feature_cat = rownames(factor_labeling)[factor_labeling$type == feature_cat]
  cat_genes_vec = as.integer(sample_set_gene_ids %in% gene_ids_in_feature_cat)
  pwf$DEgenes = cat_genes_vec
  res = goseq(pwf,gene2cat=GO_info_listed, use_genes_without_cat=FALSE)
  ## over-represented categories:
  pvals = res$over_represented_pvalue
  pvals[pvals > 1 - 1e-10] = 1 - 1e-10
  q = qvalue(pvals)
  res$over_represented_FDR = q$qvalues
  go_enrich_filename = paste(feature_cat,'.GOseq.enriched', sep='')
  result_table = res[res$over_represented_pvalue<=0.05,]
  descr = unlist(lapply(result_table$category, get_GO_term_descr))
  result_table$go_term = descr;
  result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) { 
    gene_list = GO_to_gene_list[[x]]
    gene_list = gene_list[gene_list %in% gene_ids_in_feature_cat]
    paste(gene_list, collapse=', ');
  }) )
go_list[[feature_cat]] <- result_table[order(result_table$under_represented_pvalue),]
}
View(go_list[[5]])
head(go_list[["c5"]])
