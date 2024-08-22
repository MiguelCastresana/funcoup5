# Install necessary packages if not already installed
install_if_missing <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) install.packages(package_name)
}

install_if_missing("BiocManager")
BiocManager::install("supraHex")
BiocManager::install("Rgraphviz")
install_if_missing("dnet")
install_if_missing("PRROC")

library(supraHex)
library(graph)
library(Rgraphviz)
library(dnet)
library(PRROC)
library(biomaRt)
library(igraph)
library(dplyr)

load_required_packages <- function(packages) {
  for (pkg in packages) {
    install_if_missing(pkg)
    library(pkg, character.only = TRUE)
  }
}

# create a randomized graph
create_randomized_graph <- function(df, cols = 1:2, directed = FALSE, niter_multiplier = 100) {
  df[, cols] <- lapply(df[, cols], as.vector)
  graph <- graph_from_data_frame(df[, cols], directed = directed)
  rewire(graph, with = keeping_degseq(niter = nrow(df) * niter_multiplier))
}

#  prepare seed data
prepare_seed_data <- function(genes, seeds) {
  pos <- which(genes %in% seeds)
  seed_vector <- rep(0, length(genes))
  seed_vector[pos] <- 1
  seed_data <- as.data.frame(seed_vector)
  rownames(seed_data) <- genes
  seed_data
}

# calculate the affinity matrix using dRWR
calculate_affinity_matrix <- function(graph, seed_data, restart_prob = 0.75, normalize_method = "laplacian") {
  dRWR(g = graph, normalise = normalize_method, setSeeds = seed_data, restart = restart_prob, parallel = TRUE)
}

# generate precision-recall and ROC curves
generate_curves <- function(PTmatrix, genes, end_points) {
  values <- as.data.frame(PTmatrix[, 1])
  values$Gene <- genes
  newdata <- values[order(-values[, 1]), ]
  
  true <- newdata$V1[which(newdata$Gene %in% end_points)]
  false <- newdata$V1[which(!newdata$Gene %in% end_points)]
  
  pr <- pr.curve(scores.class0 = true, scores.class1 = false)
  rc <- roc.curve(scores.class0 = true, scores.class1 = false, curve = TRUE)
  
  list(roc_auc = rc$auc, pr_auc = pr$auc.integral)
}

# perform all splits and calculate curves
all_splits <- function(KEGGA, KEGGB, funcoup_graph, genes) {
  seeds <- as.vector(KEGGA[, 1])
  end_points <- as.vector(KEGGB[, 1])
  
  seed_data <- prepare_seed_data(genes, seeds)
  PTmatrix <- calculate_affinity_matrix(funcoup_graph, seed_data)
  
  generate_curves(PTmatrix, genes, end_points)
}

# load Ensembl data and get mapping
load_ensembl_data <- function(ensembl_dataset = "hsapiens_gene_ensembl") {
  ensembl <- useMart("ensembl")
  ensembl <- useDataset(ensembl_dataset, mart = ensembl)
  getBM(attributes = c('entrezgene_id', 'ensembl_gene_id'), mart = ensembl)
}

# loop for network analysis
perform_network_analysis <- function(funcoup_f, data_createdA, data_createdB, num_networks = 30, num_splits = 30, save_path = "/scratch/Funcoup5/benchmark/results_random_networks_funcoup/") {
  roc_funcoup_true_list <- vector("list", num_splits)
  funcoup_pr_true_list <- vector("list", num_splits)
  
  for (network_index in seq_len(num_networks)) {
    funcoup_graph <- create_randomized_graph(funcoup_f)
    genes <- unique(c(as.vector(V(funcoup_graph)$name)))
    
    for (split_index in seq_len(num_splits)) {
      partA <- split(data_createdA[[split_index]], f = data_createdA[[split_index]][, 2])
      partB <- split(data_createdB[[split_index]], f = data_createdB[[split_index]][, 2])
      
      res <- mapply(all_splits, partA, partB, MoreArgs = list(funcoup_graph = funcoup_graph, genes = genes))
      
      roc_funcoup_true_list[[split_index]] <- res["roc_auc", ]
      funcoup_pr_true_list[[split_index]] <- res["pr_auc", ]
    }
    
    save(funcoup_pr_true_list, file = paste0(save_path, "funcoup_pr_true_list_network_", network_index))
    save(roc_funcoup_true_list, file = paste0(save_path, "roc_funcoup_true_list_network_", network_index))
  }
  
  list(mean_roc = mean(unlist(roc_funcoup_true_list)), mean_pr = mean(unlist(funcoup_pr_true_list)))
}


# Load required packages
required_packages <- c("supraHex", "graph", "Rgraphviz", "dnet", "PRROC", "biomaRt", "igraph", "dplyr")
load_required_packages(required_packages)

# Perform network analysis
results <- perform_network_analysis(funcoup_f, data_createdA, data_createdB, num_networks = 30, num_splits = 30)
print(paste("Mean ROC for funcoup:", results$mean_roc))
print(paste("Mean PR for funcoup:", results$mean_pr))
