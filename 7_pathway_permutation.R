#------------------------------------------------------------------------------#
#                                                                              #
#                                                                              #
#                                                                              #
#                SCRIPT TO COLLECT AND FILTER FUMA RESULTS                     #
#               SECOND PART IS FOR CYTOSCAPE NETWORK ANALYSIS                  #
#                                                                              #
#                                                                              #
#------------------------------------------------------------------------------#

# ============================================
# SETUP
# ============================================

library(dplyr)
library(tidyr)
library(biomaRt)
library(rWikiPathways)
library(RCy3)
library(ggplot2)

jaccard_similarity <- function(A, B) {
  intersection = length(intersect(A, B))
  union = length(A) + length(B) - intersection
  return (intersection/union)
}

DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(DATA.DIR)

# ============================================
# MAP TO ENTREZ GENE IDENTIFIERS
# ============================================

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- biomaRt::getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
                          mart = mart)

mapping <- mapping[!is.na(mapping$entrezgene_id),]

# ============================================
# LOAD PATHWAY INFO
# ============================================

kegg <- rWikiPathways::readGMT(file.path(getwd(), "databases", "c2.cp.kegg_medicus.v2023.2.Hs.entrez.gmt"))
reactome <- rWikiPathways::readGMT(file.path(getwd(), "databases", "c2.cp.reactome.v2023.2.Hs.entrez.gmt"))
wp <- rWikiPathways::readGMT(file.path(getwd(), "databases", "c2.cp.wikipathways.v2023.2.Hs.entrez.gmt"))
pathways <- rbind(kegg, reactome, wp)
rm(kegg, reactome, wp)

pathways <- pathways[pathways$gene %in% mapping$entrezgene_id,]

mapping <- mapping[mapping$entrezgene_id %in% pathways$gene,]

# ============================================
# RANDOM SAMPLING OVERLAP IN GENES & PATHWAYS
# ============================================

# JACCARD INDEX OF GENES IN PATHWAYS
# COLLECTED FROM RANDOM GENE SETS
# ALSO GET THE TOTAL GENE SET SIZES
# AND THE TOTAL OF ASSOCIATED PATHWAYS
overlap_genes <- overlap_paths <- data.frame(matrix(nrow = 100000, ncol = 3))
colnames(overlap_genes) <- colnames(overlap_paths) <- c("pw_ps", "pw_ls", "ps_ls")

total_pathways <- total_set_size <- data.frame(matrix(nrow = 100000, ncol = 3))
colnames(total_pathways) <- colnames(total_set_size) <- c("pw", "ps", "ls")

unique_total_size <- data.frame(matrix(nrow = 100000, ncol = 3))
colnames(unique_total_size) <- c("pw_ps", "pw_ls", "ps_ls")

set.seed(14021996)

for (i in 1:100000) {
  # Random sampling from genes in pathways
  pw <- mapping[sample(nrow(mapping), 19),]
  ps <- mapping[sample(nrow(mapping), 22),]
  ls <- mapping[sample(nrow(mapping), 8),]
  
  # Take those sets
  pw <- pathways[pathways$gene %in% pw$entrezgene_id,]
  ps <- pathways[pathways$gene %in% ps$entrezgene_id,]
  ls <- pathways[pathways$gene %in% ls$entrezgene_id,]
  
  # Get full pathway gene sets from those selected pathways
  pw_asset <- pathways[pathways$term %in% pw$term,]
  ps_asset <- pathways[pathways$term %in% ps$term,]
  ls_asset <- pathways[pathways$term %in% ls$term,]
  
  # Get total associated gene set size
  # as selected by the random gene set pathways
  total_set_size[i,1] <- length(unique(pw_asset$gene))
  total_set_size[i,2] <- length(unique(ps_asset$gene))
  total_set_size[i,3] <- length(unique(ls_asset$gene))
  
  # Get total associated pathways 
  # as selected by the random gene set 
  total_pathways[i,1] <- length(unique(pw$term))
  total_pathways[i,2] <- length(unique(ps$term))
  total_pathways[i,3] <- length(unique(ls$term))
  
  # Get unique associated genes
  pw_genes <- unique(pw$gene)
  ps_genes <- unique(ps$gene)
  ls_genes <- unique(ls$gene)
  
  # Get total unique associated gene set size
  unique_total_size[i,1] <- length(unique(c(unique(pw_asset$gene), unique(ps_asset$gene))))
  unique_total_size[i,2] <- length(unique(c(unique(pw_asset$gene), unique(ls_asset$gene))))
  unique_total_size[i,3] <- length(unique(c(unique(ps_asset$gene), unique(ls_asset$gene))))
  
  # Get Jaccard similarity of the associated gene sets
  overlap_genes[i,1] <- jaccard_similarity(pw_asset$gene, ps_asset$gene)
  overlap_genes[i,2] <- jaccard_similarity(pw_asset$gene, ls_asset$gene)
  overlap_genes[i,3] <- jaccard_similarity(ps_asset$gene, ls_asset$gene)
  
  # Get unique associated pathways
  pw_paths <- unique(pw$term)
  ps_paths <- unique(ps$term)
  ls_paths <- unique(ls$term)
  
  # Get Jaccard similarity of the pathways
  overlap_paths[i,1] <- jaccard_similarity(pw_paths, ps_paths)
  overlap_paths[i,2] <- jaccard_similarity(pw_paths, ls_paths)
  overlap_paths[i,3] <- jaccard_similarity(ps_paths, ls_paths)
}

hist(overlap_paths$pw_ps/unique_total_size[,1], breaks = 250,
     xlab = "Pathway overlap / by total gene set size", 
     ylab = "Frequency", xlim = c(0,7e-5), ylim = c(0,2500),
     main = "F6 A - Penn Word with List Sorting")
abline(v=1.304123e-05,col="red")
# percentile 44

hist(overlap_paths$pw_ls/unique_total_size[,2], breaks = 500,
     xlab = "Pathway overlap / by total gene set size", 
     ylab = "Frequency", xlim = c(0,9e-5), ylim = c(0,3000),
     main = "F6 B - Penn Word with Picture Sequence")
abline(v=1.422275e-05,col="red")
# percentile 76

hist(overlap_paths$ps_ls/unique_total_size[,3], breaks = 500,
     xlab = "Pathway overlap / by total gene set size", 
     ylab = "Frequency", xlim = c(0,8e-5), ylim = c(0,4500),
     main = "F6 C - Picture Sequence with List Sorting")
abline(v=1.420937e-05,col="red")
# percentile 80

quantile(overlap_paths$pw_ps/unique_total_size[,1], probs = seq(0, 1, 0.01))

quantile(na.omit(overlap_paths$pw_ls/unique_total_size[,2]), probs = seq(0, 1, 0.01))

quantile(overlap_paths$ps_ls/unique_total_size[,3], probs = seq(0, 1, 0.01))





