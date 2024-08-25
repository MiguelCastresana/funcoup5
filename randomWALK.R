# Load necessary libraries
library(tidyverse)
library(mygene)
library(igraph)
library(dnet)
library(PRROC)

# Define custom negation function
`%!in%` <- Negate(`%in%`)

# Load and process the 'string' dataset
string <- read_delim("~/FC5.0_H.sapiens_compact", delim = "\t") %>%
  filter(X1 >= 0.7) %>%
  select(X3, X4, X1)

# Create degree list
degreelist <- string %>%
  select(X3, X4) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  count(gene, name = "degree") %>%
  arrange(desc(degree))

# Extract the genome as a vector
genome <- degreelist %>% pull(gene)

# Load input data
input_data <- read_delim("~/orphanet_greater10", delim = "\t")

# Extract unique genes
genes <- input_data %>% pull(X1) %>% unique()

# Query gene information
out <- queryMany(genes, scopes = "symbol", fields = "ensembl.gene", species = "human")

# Process the output data
dat_f <- map_df(seq_along(out$query), function(i) {
  ok <- out$ensembl[[i]]$gene
  tibble(values = rep(out$query[[i]], length(ok)), ensembl = ok)
})

# Join data and filter
final <- dat_f %>%
  rename(group = ensembl) %>%
  inner_join(input_data, by = c("values" = "X1")) %>%
  select(ensembl, group) %>%
  distinct() %>%
  filter(ensembl %in% genome)

# Further filtering based on selected groups
filtered_data <- final %>% filter(group %in% selected)

iterations = 30

# Prepare datasets for analysis
data_created <- map(1:iterations, function(times) {
  groups <- unique(filtered_data$group)
  final_dat1 <- tibble()
  final_dat2 <- tibble()
  
  for (grp in groups) {
    genes <- filtered_data %>%
      filter(group == grp) %>%
      pull(ensembl)
    
    if (length(genes) < 5) next
    
    g1 <- round(length(genes) / 2)
    sub1 <- sample(genes, g1)
    sub2 <- setdiff(genes, sub1)
    
    final_dat1 <- bind_rows(final_dat1, tibble(sub1, group = grp))
    final_dat2 <- bind_rows(final_dat2, tibble(sub2, group = grp))
  }
  
  list(data_A = final_dat1, data_B = final_dat2)
})

# Save datasets
save(map(data_created, "data_A"), file = "/scratch/Funcoup5/orphanet/data/splitA_personalized_07")
save(map(data_created, "data_B"), file = "/scratch/Funcoup5/orphanet/data/splitB_personalized_07")

# Function for performing splits and calculating PR and ROC curves
all_splits <- function(part1, part2) {
  seeds <- part1$ensembl
  end_points <- part2$ensembl
  
  seed_vector <- as.integer(genes_graph %in% seeds)
  seed_data <- as_tibble(seed_vector, .name_repair = "minimal")
  rownames(seed_data) <- genes_graph
  
  PTmatrix <- dRWR(g = graph_chosen, normalise = "laplacian", setSeeds = seed_data, restart = 0.5, parallel = TRUE)
  
  values <- as_tibble(PTmatrix[,1], .name_repair = "minimal") %>%
    bind_cols(genes_graph) %>%
    arrange(desc(value))
  
  true_scores <- values %>% filter(value %in% end_points) %>% pull(value)
  false_scores <- values %>% filter(value %!in% end_points) %>% pull(value)
  
  pr <- pr.curve(scores.class0 = true_scores, scores.class1 = false_scores)
  rc <- roc.curve(scores.class0 = true_scores, scores.class1 = false_scores)
  
  return(c(rc$auc, pr$auc.integral))
}

# Build the string graph
string_graph <- graph_from_data_frame(d = string[, 1:2], directed = FALSE)
genes_graph <- V(string_graph)$name

# Calculate true PR and ROC curves
results <- map(1:iterations, function(splits) {
  partA <- split(data_created[[splits]]$data_A, f = data_created[[splits]]$data_A$group)
  partB <- split(data_created[[splits]]$data_B, f = data_created[[splits]]$data_B$group)
  
  res <- mapply(all_splits, partA, partB, SIMPLIFY = FALSE)
  
  list(roc = res[seq(1, length(res), 2)], pr = res[seq(2, length(res), 2)])
})

# Save final results
save(map(results, "pr"), file = "/scratch/Funcoup5/orphanet/results/TRUE_pr_ALL_random_walk_personalized_07")
save(map(results, "roc"), file = "/scratch/Funcoup5/orphanet/results/TRUE_roc_ALL_random_walk_personalized_07")
