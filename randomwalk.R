
# Random walk with restart through three different networks: FunCoup, STRING and HumanNet

# Load necessary libraries
library(dnet)
library(PRROC)
library(igraph)
library(biomaRt)

# Optional: install packages only if not installed
# install.packages("PRROC")
# BiocManager::install("dnet")

# Load network data (e.g., FunCoup, STRING, HumanNet)
load_network <- function(path, sep = "\t", header = TRUE) {
  network <- read.delim(path, sep = sep, header = header)
  return(network[, 1:2])  # Ensure only edge pairs are returned
}

# Function to calculate PR and ROC
run_affinity <- function(graph, seeds, end_points, genes) {
  pos <- which(genes %in% seeds)
  pos_end <- which(genes %in% end_points)
  
  seed_vector <- rep(0, length(genes))
  seed_vector[pos] <- 1
  seed_data <- as.data.frame(seed_vector)
  rownames(seed_data) <- genes
  
  PTmatrix <- dRWR(g = graph, normalise = "laplacian", setSeeds = seed_data, restart = 0.75, parallel = TRUE)
  values <- as.data.frame(PTmatrix[, 1])
  values[, 2] <- genes
  newdata <- values[order(-values[, 1]), ]
  
  true <- newdata[which(newdata[, 2] %in% end_points), 1]
  false <- newdata[which(newdata[, 2] %!in% end_points), 1]
  
  pr <- pr.curve(scores.class0 = true, scores.class1 = false)
  rc <- roc.curve(scores.class0 = true, scores.class1 = false)
  
  return(list(pr_auc = pr$auc.integral, roc_auc = rc$auc))
}

# Create KEGGA/KEGGB splits
create_splits <- function(KEGG, n = 30) {
  groups <- unique(as.vector(KEGG[, 2]))
  dataA <- list()
  dataB <- list()
  
  for (split in 1:n) {
    final_dat1 <- data.frame()
    final_dat2 <- data.frame()
    
    for (grp in groups) {
      genes <- as.vector(KEGG[KEGG[, 2] == grp, 1])
      g1 <- sample(genes, round(length(genes) / 2))
      g2 <- setdiff(genes, g1)
      
      dat1 <- data.frame(ensembl_gene_id = g1, disease = grp)
      dat2 <- data.frame(ensembl_gene_id = g2, disease = grp)
      
      final_dat1 <- rbind(final_dat1, dat1)
      final_dat2 <- rbind(final_dat2, dat2)
    }
    
    dataA[[split]] <- final_dat1
    dataB[[split]] <- final_dat2
  }
  
  return(list(A = dataA, B = dataB))
}

# Benchmark network with KEGG-based seed splits
benchmark_network <- function(graph_data, KEGG, splits, output_prefix) {
  results_pr <- list()
  results_roc <- list()
  
  for (network_id in 1:30) {
    graph <- graph_from_data_frame(graph_data[, 1:2], directed = FALSE)
    graph <- rewire(graph, with = keeping_degseq(niter = nrow(graph_data) * 100))
    genes <- unique(c(as.vector(graph_data[, 1]), as.vector(graph_data[, 2])))
    V(graph)$name <- genes
    
    pr_list <- numeric()
    roc_list <- numeric()
    
    for (i in 1:length(splits$A)) {
      partA <- split(splits$A[[i]], splits$A[[i]][, 2])
      partB <- split(splits$B[[i]], splits$B[[i]][, 2])
      
      for (group in names(partA)) {
        res <- run_affinity(graph, partA[[group]][, 1], partB[[group]][, 1], genes)
        pr_list <- c(pr_list, res$pr_auc)
        roc_list <- c(roc_list, res$roc_auc)
      }
    }
    
    results_pr[[network_id]] <- pr_list
    results_roc[[network_id]] <- roc_list
    
    save(results_pr[[network_id]], file = paste0(output_prefix, "_pr_network_", network_id))
    save(results_roc[[network_id]], file = paste0(output_prefix, "_roc_network_", network_id))
  }
}

# Example usage:
# KEGG <- read.delim("/path/to/your/kegg_file.tsv")
# splits <- create_splits(KEGG)
# funcoup_f <- load_network("/path/to/funcoup_network.tsv")
# benchmark_network(funcoup_f, KEGG, splits, "/output/funcoup")
