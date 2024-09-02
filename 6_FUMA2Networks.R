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
library(biomaRt)
if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways")
}
library(rWikiPathways)
library(RCy3)
library(tidyr)
library(tibble)

DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(DATA.DIR)

extract_shared_SNP_interactions <- function(x) {
  result <- x %>% group_by(uniqID) %>% filter(n() > 1) %>%
    summarise(gene = toString(unique(gene))) %>% ungroup()
  
  
  generate_pairs <- function(x) {
    elements <- strsplit(x, ", ")[[1]]
    combn(elements, 2, simplify = FALSE)
  }
  
  # Generate the pairs and create a new data frame
  pairs_df <- unique(result %>% rowwise() %>%
                       mutate(pairs = list(generate_pairs(gene))) %>%
                       unnest(cols = c(pairs)) %>%
                       mutate(source = sapply(pairs, `[`, 1),
                              target = sapply(pairs, `[`, 2)) %>%
                       dplyr::select(source, target))
  
  return(pairs_df)
}

# ============================================
# READ FUMA OUTPUT
# ============================================

# Penn word 
# ============================================
# get gene information of SNPs located in the genes
pw.annov <- read.table(paste0("FUMA_iwrd/annov.txt"), header = T, sep = "\t")
pw.annov.filt <- pw.annov[pw.annov$annot %in% c("intronic","exonic","UTR3","UTR5"),]
pw.annov.genes <- pw.annov.filt[,c("gene", "symbol")]
pw.annov.genes$memory.test <- "PW"
pw.annov.genes$type <- "SNP"
pw.annov.genes <- unique(pw.annov.genes)

pw.annov.snp.gene <- pw.annov.filt[,c(1,2)]
pw.annov.snp.gene <- unique(pw.annov.snp.gene)

rm(pw.annov, pw.annov.filt)

# get gene information of eQTLs
pw.eqtl <- read.table(paste0("FUMA_iwrd/eqtl.txt"), header = T, sep = "\t")
pw.eqtl.genes <- pw.eqtl[,c("gene", "symbol")]
pw.eqtl.genes$memory.test <- "PW"
pw.eqtl.genes$type <- "eQTL"
pw.eqtl.genes <- unique(pw.eqtl.genes)

# get information about shared SNPs 
pw.eqtl.snp.gene <- pw.eqtl[,c(1,4)]
pw.snp.gene <- rbind(pw.annov.snp.gene, pw.eqtl.snp.gene)
pw.snp.gene <- unique(pw.snp.gene)

pw.genes <- rbind(pw.annov.genes, pw.eqtl.genes)
pw.genes <- pw.genes %>% group_by(gene, symbol, memory.test) %>% summarise(type = toString(unique(type))) 
pw.shared.snps <- extract_shared_SNP_interactions(pw.snp.gene)


rm(pw.eqtl, pw.eqtl.snp.gene, pw.annov.snp.gene, pw.snp.gene, pw.annov.genes, pw.eqtl.genes)

# List sorting
# ============================================
# get gene information of SNPs located in the genes
ls.annov <- read.table(paste0("FUMA_list/annov.txt"), header = T, sep = "\t")
ls.annov.filt <- ls.annov[ls.annov$annot %in% c("intronic","exonic","UTR3","UTR5"),]
ls.annov.genes <- ls.annov.filt[,c("gene", "symbol")]
ls.annov.genes$memory.test <- "LS"
ls.annov.genes$type <- "SNP"
ls.annov.genes <- unique(ls.annov.genes)

ls.annov.snp.gene <- ls.annov.filt[,c(1,2)]
ls.annov.snp.gene <- unique(ls.annov.snp.gene)

rm(ls.annov, ls.annov.filt)

# get gene information of eQTLs
ls.eqtl <- read.table(paste0("FUMA_list/eqtl.txt"), header = T, sep = "\t")
ls.eqtl.genes <- ls.eqtl[,c("gene", "symbol")]
ls.eqtl.genes$memory.test <- "LS"
ls.eqtl.genes$type <- "eQTL"
ls.eqtl.genes <- unique(ls.eqtl.genes)

# get information about shared SNPs 
ls.eqtl.snp.gene <- ls.eqtl[,c(1,4)]
ls.snp.gene <- rbind(ls.annov.snp.gene, ls.eqtl.snp.gene)
ls.snp.gene <- unique(ls.snp.gene)

ls.genes <- rbind(ls.annov.genes, ls.eqtl.genes)
ls.genes <- ls.genes %>% group_by(gene, symbol, memory.test) %>% summarise(type = toString(unique(type))) 
ls.shared.snps <- extract_shared_SNP_interactions(ls.snp.gene)


rm(ls.eqtl, ls.eqtl.snp.gene, ls.annov.snp.gene, ls.snp.gene, ls.annov.genes, ls.eqtl.genes)


# Picture sequence
# ============================================
# get gene information of SNPs located in the genes
ps.annov <- read.table(paste0("FUMA_pics/annov.txt"), header = T, sep = "\t")
ps.annov.filt <- ps.annov[ps.annov$annot %in% c("intronic","exonic","UTR3","UTR5"),]
ps.annov.genes <- ps.annov.filt[,c("gene", "symbol")]
ps.annov.genes$memory.test <- "PS"
ps.annov.genes$type <- "SNP"
ps.annov.genes <- unique(ps.annov.genes)

ps.annov.snp.gene <- ps.annov.filt[,c(1,2)]
ps.annov.snp.gene <- unique(ps.annov.snp.gene)

rm(ps.annov, ps.annov.filt)

# get gene information of eQTps
ps.eqtl <- read.table(paste0("FUMA_pics/eqtl.txt"), header = T, sep = "\t")
ps.eqtl.genes <- ps.eqtl[,c("gene", "symbol")]
ps.eqtl.genes$memory.test <- "PS"
ps.eqtl.genes$type <- "eQTL"
ps.eqtl.genes <- unique(ps.eqtl.genes)

# get information about shared SNPs 
ps.eqtl.snp.gene <- ps.eqtl[,c(1,4)]
ps.snp.gene <- rbind(ps.annov.snp.gene, ps.eqtl.snp.gene)
ps.snp.gene <- unique(ps.snp.gene)

ps.genes <- rbind(ps.annov.genes, ps.eqtl.genes)
ps.genes <- ps.genes %>% group_by(gene, symbol, memory.test) %>% summarise(type = toString(unique(type))) 
ps.shared.snps <- extract_shared_SNP_interactions(ps.snp.gene)


rm(ps.eqtl, ps.eqtl.snp.gene, ps.annov.snp.gene, ps.snp.gene, ps.annov.genes, ps.eqtl.genes)


# ============================================
# COMBINE DIFFERENT MEMORY TESTS
# ============================================

# MERGE AND ONLY KEEP UNIQUE GENES
genes <- rbind(pw.genes, ls.genes, ps.genes)
shared.snps <- rbind(pw.shared.snps,ls.shared.snps, ps.shared.snps)
shared.snps$interaction <- "shared SNP"
# ============================================
# MAP TO ENTREZ GENE IDENTIFIERS
# ============================================

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- biomaRt::getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","entrezgene_id"),
               values = unique(genes$gene), 
               mart = mart)

genes <- merge(genes, mapping, by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE)
genes$node.type <- "Gene"
genes <- genes[order(genes$entrezgene_id), ]
genes <- genes[!duplicated(genes$gene), ]
genes$entrezgene_id[is.na(genes$entrezgene_id)] <- genes$symbol[is.na(genes$entrezgene_id)]

shared.snps$source <- genes$entrezgene_id[match(shared.snps$source, genes$gene)]
shared.snps$target <- genes$entrezgene_id[match(shared.snps$target, genes$gene)]
shared.snps <- as.data.frame(t(apply(shared.snps, 1, function(x) sort(x))))
shared.snps <- unique(shared.snps)
colnames(shared.snps) <- c("source", "target", "interaction")

# ============================================
# LOAD PATHWAY INFO
# ============================================

kegg <- rWikiPathways::readGMT(file.path(getwd(), "databases", "c2.cp.kegg_medicus.v2023.2.Hs.entrez.gmt"))
reactome <- rWikiPathways::readGMT(file.path(getwd(), "databases", "c2.cp.reactome.v2023.2.Hs.entrez.gmt"))
wp <- rWikiPathways::readGMT(file.path(getwd(), "databases", "c2.cp.wikipathways.v2023.2.Hs.entrez.gmt"))
pathways <- rbind(kegg, reactome, wp)
rm(kegg, reactome, wp)

pathways.filt <- pathways[pathways$gene %in% genes$entrezgene_id,]
colnames(pathways.filt) <- c("source", "target")
pathways.filt$interaction <- "pathway-gene"

nodes.g <- genes[,c(5,2)]
colnames(nodes.g) <- c("id", "label")
nodes.p <- unique(pathways.filt[,c(1,1)])
colnames(nodes.p) <- c("id", "label")
nodes <- rbind(nodes.g, nodes.p)

nodes.p$node.type <- "pathway"

edges <- rbind(pathways.filt, shared.snps)

# ============================================
# CREATE CYTOSCAPE NETWORK
# ============================================

RCy3::createNetworkFromDataFrames(nodes=nodes,edges=edges, title="FUMA-all-memory-tests-merged", collection="HPA-FUMA")
RCy3::loadTableData(genes, data.key.column = "entrezgene_id", table.key.column = "id")
RCy3::loadTableData(nodes.p, data.key.column = "label", table.key.column = "label")

RCy3::createVisualStyle("pathway-gene")
RCy3::setNodeColorMapping(table.column = "memory.test",
                    table.column.values = c("PW", "LS", "PS"),
                    mapping.type = "d",
                    colors = c("#ff9933", "#6699cc", "#99cc99"), default.color = "#D9D9D9", style.name = "pathway-gene")
RCy3::setNodeShapeMapping("node.type", table.column.values = c("Gene", "pathway"), shapes = c("ELLIPSE","DIAMOND"), style.name = "pathway-gene")
RCy3::setNodeSizeMapping("node.type", table.column.values = c("Gene", "pathway"), sizes = c(60,40), mapping.type = "d",style.name = "pathway-gene")
RCy3::setNodeLabelMapping("symbol", style.name = "pathway-gene")
RCy3::setNodeFontSizeDefault("25", style.name = "pathway-gene")
RCy3::setNodeBorderWidthDefault(0,style.name = "pathway-gene")
RCy3::setEdgeColorDefault("#AFAAAA",style.name = "pathway-gene")
RCy3::setEdgeLineWidthDefault(2,style.name = "pathway-gene")
RCy3::setEdgeLineStyleMapping("interaction", c("pathway-gene", "shared SNP"), c("SOLID", "LONG_DASH"), "SOLID", "pathway-gene" )

RCy3::setVisualStyle("pathway-gene")

RCy3::toggleGraphicsDetails()

# Manual layouting!

# ============================================
# FILTER NETWORK BASED ON PATHWAYS WITH MULTIPLE GENES
# ============================================

toCamelCase <- function(s) {
  allCapsWords <- c("RNA", "DNA", "TP53", "NRF2", "GPCR", "GTPase", "GTPases", "SLC", "RHOBTB3")
  # Split the string into words based on underscores or spaces
  words <- unlist(strsplit(s, "_| "))
  
  words <- sapply(words, function(word) {
    # If the word is in the list of allCapsWords, leave it in all caps
    if (toupper(word) %in% toupper(allCapsWords)) {
      return(toupper(word))
    } else {
      # Otherwise, convert it to CamelCase
      return(paste0(toupper(substring(word, 1, 1)), tolower(substring(word, 2))))
    }
  })
  
  # Combine the words into a single string
  camelCaseString <- paste0(words, collapse = " ")
  return(camelCaseString)
}

RCy3::analyzeNetwork(directed=FALSE)
df.degree <- RCy3::getTableColumns("node", c("SUID","id","Degree", "node.type"))
df.degree.pwys <- subset(df.degree, node.type=="pathway" & Degree >= 2)
selected.pwys <- df.degree.pwys$id
df.degree.pwys <- df.degree.pwys %>% separate(id, into = c("type", "label"), sep = "_", extra = "merge")
df.degree.pwys$label <- sapply(df.degree.pwys$label,toCamelCase)
df.degree.genes <- subset(df.degree, node.type=="Gene")

RCy3::loadTableData(data = df.degree.pwys, data.key.column = "SUID", table = "node", table.key.column ="SUID")
select.nodes <- rbind(df.degree.pwys[,c("SUID","node.type")], df.degree.genes[,c("SUID","node.type")])
RCy3::createSubnetwork(nodes=select.nodes$SUID, subnetwork.name = "all tests pathways > 2 genes")
RCy3::layoutNetwork()
RCy3::cloneNetwork()

df.edges <- RCy3::getTableColumns("edge", c("SUID","interaction"))
df.edges <- subset(df.edges, interaction=="shared SNP")
RCy3::selectEdges(df.edges$SUID)
RCy3::deleteSelectedEdges()

RCy3::createVisualStyle("overlap")
RCy3::setNodeColorMapping(table.column = "memory.test",
                          table.column.values = c("PW", "LS", "PS"),
                          mapping.type = "d",
                          colors = c("#ff9933", "#6699cc", "#99cc99"), default.color = "#D9D9D9", style.name = "overlap")
RCy3::setNodeShapeMapping("node.type", table.column.values = c("Gene", "pathway"), shapes = c("ELLIPSE","DIAMOND"), style.name = "overlap")
RCy3::setNodeSizeMapping("node.type", table.column.values = c("Gene", "pathway"), sizes = c(50,25), mapping.type = "d",style.name = "overlap")
RCy3::setNodeLabelMapping("label", style.name = "overlap")
RCy3::setNodeFontSizeMapping("node.type", table.column.values = c("Gene", "pathway"), sizes = c(12,20), mapping.type = "d",style.name = "overlap")
RCy3::setNodeFontSizeDefault("25", style.name = "overlap")
RCy3::setNodeBorderWidthDefault(0,style.name = "overlap")
RCy3::setEdgeColorDefault("#AFAAAA",style.name = "overlap")
RCy3::setEdgeLineWidthDefault(2,style.name = "overlap")
RCy3::setEdgeLineStyleMapping("interaction", c("pathway-gene", "shared SNP"), c("SOLID", "LONG_DASH"), "SOLID", "overlap" )

RCy3::setVisualStyle("overlap")

RCy3::toggleGraphicsDetails()

pathway_overlap <- function(pathways, selected.pathways) {
  overlap <- data.frame(matrix(ncol = 5))
  colnames(overlap) <- c("source", "target", "size.source", "size.target", "max.overlap")
  
  n.pwys <- length(selected.pathways)
  
  x <- 1
  for(i in 1:n.pwys) {
    for(j in (i+1):n.pwys) {
      if(i != j) {
        genes.x <- pathways[pathways$term==selected.pathways[i],2]
        genes.y <- pathways[pathways$term==selected.pathways[j],2]
        
        is <- intersect(genes.x, genes.y)
        o1 <- length(is)/length(genes.x)
        o2 <- length(is)/length(genes.y)
        if(o1 > 0.8 || o2 > 0.8) {
          overlap[x,1] <- selected.pathways[i]
          overlap[x,2] <- selected.pathways[j]
          overlap[x,3] <- length(genes.x)
          overlap[x,4] <- length(genes.y)
          overlap[x,5] <- max(o1,o2)
          
          x <- x + 1
        } 
      }
    }
  }
  return(overlap)
}

overlap <- pathway_overlap(pathways,selected.pwys)

## MANUALLY CHECK OVERLAP AND SEE WHICH PATHWAYS SHOULD BE MERGED
nodes.cellcycle <- c("WP_CELL_CYCLE", "REACTOME_CELL_CYCLE","REACTOME_CELL_CYCLE_MITOTIC","REACTOME_CELL_CYCLE_CHECKPOINTS")
RCy3::createGroup("Cell cycle", nodes = nodes.cellcycle, nodes.by.col="id")
RCy3::collapseGroup("Cell cycle")
df.cellcycle <- tibble(id = "Cell cycle", label = "Cell cycle", type = "group", node.type = "pathway", pathways = paste(nodes.cellcycle, collapse=","))
RCy3::loadTableData(df.cellcycle, data.key.column = "id", table = "node", table.key.column = "shared name")


nodes.apoptosis <- c("REACTOME_ACTIVATION_OF_BH3_ONLY_PROTEINS", "REACTOME_APOPTOSIS","REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS","REACTOME_ACTIVATION_OF_PUMA_AND_TRANSLOCATION_TO_MITOCHONDRIA", "REACTOME_PROGRAMMED_CELL_DEATH", "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53")
RCy3::createGroup("Apoptosis", nodes = nodes.apoptosis, nodes.by.col="id")
RCy3::collapseGroup("Apoptosis")
df.apoptosis <- tibble(id = "Apoptosis", label = "Apoptosis", type = "group", node.type = "pathway", pathways = paste(nodes.apoptosis, collapse=","))
RCy3::loadTableData(df.apoptosis, data.key.column = "id", table = "node", table.key.column = "shared name")


nodes.nr <- c("WP_NRF2_PATHWAY", "WP_NUCLEAR_RECEPTORS_META_PATHWAY")
RCy3::createGroup("Nuclear receptors", nodes = nodes.nr, nodes.by.col="id")
RCy3::collapseGroup("Nuclear receptors")
df.nr <- tibble(id = "Nuclear receptors", label = "Nuclear receptors", type = "group", node.type = "pathway", pathways = paste(nodes.nr, collapse=","))
RCy3::loadTableData(df.nr, data.key.column = "id", table = "node", table.key.column = "shared name")

nodes.cytokine <- c("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM", "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING","REACTOME_INTERFERON_SIGNALING", "WP_IMMUNE_RESPONSE_TO_TUBERCULOSIS")
RCy3::createGroup("Cytokine signaling", nodes = nodes.cytokine, nodes.by.col="id")
RCy3::collapseGroup("Cytokine signaling")
df.cytokine <- tibble(id = "Cytokine signaling", label = "Cytokine signaling", type = "group", node.type = "pathway", pathways = paste(nodes.cytokine, collapse=","))
RCy3::loadTableData(df.cytokine, data.key.column = "id", table = "node", table.key.column = "shared name")

nodes.gpcr <- c("REACTOME_CLASS_A_1_RHODOPSIN_LIKE_RECEPTORS", "REACTOME_GPCR_LIGAND_BINDING","REACTOME_SIGNALING_BY_GPCR","REACTOME_G_ALPHA_Q_SIGNALLING_EVENTS")
RCy3::createGroup("GPCR signaling", nodes = nodes.gpcr, nodes.by.col="id")
RCy3::collapseGroup("GPCR signaling")
df.gpcr <- tibble(id = "GPCR signaling", label = "GPCR signaling", type = "group", node.type = "pathway", pathways = paste(nodes.gpcr, collapse=","))
RCy3::loadTableData(df.gpcr, data.key.column = "id", table = "node", table.key.column = "shared name")

nodes.contraction <- c("REACTOME_CARDIAC_CONDUCTION", "REACTOME_MUSCLE_CONTRACTION","REACTOME_PHASE_0_RAPID_DEPOLARISATION", "WP_MYOMETRIAL_RELAXATION_AND_CONTRACTION_PATHWAYS")
RCy3::createGroup("Contraction related processes", nodes = nodes.contraction, nodes.by.col="id")
RCy3::collapseGroup("Contraction related processes")
df.contraction <- tibble(id = "Contraction related processes", label = "Contraction related processes", type = "group", node.type = "pathway", pathways = paste(nodes.contraction, collapse=","))
RCy3::loadTableData(df.contraction, data.key.column = "id", table = "node", table.key.column = "shared name")

nodes.fa <- c("WP_PI3K_AKT_SIGNALING_PATHWAY", "WP_FOCAL_ADHESION_PI3K_AKT_MTOR_SIGNALING_PATHWAY")
RCy3::createGroup("Focal adhesion", nodes = nodes.fa, nodes.by.col="id")
RCy3::collapseGroup("Focal adhesion")
df.fa <- tibble(id = "Focal adhesion", label = "Focal adhesion", type = "group", node.type = "pathway", pathways = paste(nodes.fa, collapse=","))
RCy3::loadTableData(df.fa, data.key.column = "id", table = "node", table.key.column = "shared name")

nodes.rho <- c("REACTOME_RHO_GTPASE_EFFECTORS", "REACTOME_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3")
RCy3::createGroup("RHO GTPase signaling", nodes = nodes.rho, nodes.by.col="id")
RCy3::collapseGroup("RHO GTPase signaling")
df.rho <- tibble(id = "RHO GTPase signaling", label = "RHO GTPase signaling", type = "group", node.type = "pathway", pathways = paste(nodes.rho, collapse=","))
RCy3::loadTableData(df.rho, data.key.column = "id", table = "node", table.key.column = "shared name")

nodes.infection <- c("REACTOME_INFECTIOUS_DISEASE", "REACTOME_VIRAL_INFECTION_PATHWAYS")
RCy3::createGroup("Viral infection", nodes = nodes.infection, nodes.by.col="id")
RCy3::collapseGroup("Viral infection")
df.infection <- tibble(id = "Viral infection", label = "Viral infection", type = "group", node.type = "pathway", pathways = paste(nodes.infection, collapse=","))
RCy3::loadTableData(df.infection, data.key.column = "id", table = "node", table.key.column = "shared name")

nodes.transport <- c("REACTOME_TRANSPORT_OF_SMALL_MOLECULES", "REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT", "WP_PROXIMAL_TUBULE_TRANSPORT")
RCy3::createGroup("Transport of small molecules", nodes = nodes.transport, nodes.by.col="id")
RCy3::collapseGroup("Transport of small molecules")
df.transport <- tibble(id = "Transport of small molecules", label = "Transport of small molecules", type = "group", node.type = "pathway", pathways = paste(nodes.transport))
RCy3::loadTableData(df.transport, data.key.column = "id", table = "node", table.key.column = "shared name")

nodes.dev <- c("REACTOME_DEVELOPMENTAL_BIOLOGY", "REACTOME_NERVOUS_SYSTEM_DEVELOPMENT")
RCy3::createGroup("Nervous system development", nodes = nodes.dev, nodes.by.col="id")
RCy3::collapseGroup("Nervous system development")
df.dev <- tibble(id = "Nervous system development", label = "Nervous system development", type = "group", node.type = "pathway", pathways = paste(nodes.dev))
RCy3::loadTableData(df.dev, data.key.column = "id", table = "node", table.key.column = "shared name")

nodes.transport2 <- c("REACTOME_DISORDERS_OF_TRANSMEMBRANE_TRANSPORTERS", "REACTOME_SLC_TRANSPORTER_DISORDERS")
RCy3::createGroup("Transporter disorders", nodes = nodes.transport2, nodes.by.col="id")
RCy3::collapseGroup("Transporter disorders")
df.transport2 <- tibble(id = "Transporter disorders", label = "Transporter disorders", type = "group", node.type = "pathway", pathways = paste(nodes.transport2))
RCy3::loadTableData(df.transport2, data.key.column = "id", table = "node", table.key.column = "shared name")

nodes.transport2 <- c("REACTOME_NEURONAL_SYSTEM", "REACTOME_PROTEIN_PROTEIN_INTERACTIONS_AT_SYNAPSES")
RCy3::createGroup("Synapses and neuronal system", nodes = nodes.transport2, nodes.by.col="id")
RCy3::collapseGroup("Synapses and neuronal system")
df.transport2 <- tibble(id = "Synapses and neuronal system", label = "Synapses and neuronal system", type = "group", node.type = "pathway", pathways = paste(nodes.transport2))
RCy3::loadTableData(df.transport2, data.key.column = "id", table = "node", table.key.column = "shared name")

RCy3::layoutNetwork()

# FIX LAYOUT AND RENAME SOME OF THE PATHWAYS IF NEEDED