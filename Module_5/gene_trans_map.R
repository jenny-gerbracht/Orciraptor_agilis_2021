###################################################################
#Script Name	:gene_trans_map.R		                                                                                              
#Description	:Generates gene trans map for rnaSPAdes assembly                                                       
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

library(seqinr)
library(stringr)

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
