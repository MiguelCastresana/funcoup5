# Load required libraries
library(biomaRt)
library(dplyr)

# Load gold standards
ppi_FC4 <- read.delim("~/Desktop/Gold_standards/ppi_fc4")
ppi_FC5 <- read.delim("~/Desktop/Gold_standards/ppi_FC5")
ids_fc4 <- unique(ppi_FC4[, 1])
ids_fc5 <- unique(ppi_FC5[, 1])
common_ids <- intersect(ids_fc4, ids_fc5)

# Load complex datasets
complex_FC4 <- read.delim("~/Desktop/Gold_standards/complex_fc4.txt", header = FALSE, sep = ",")
complex_FC5 <- read.delim("~/Desktop/Gold_standards/complex_FC5", header = FALSE, sep = ",")

# Load and normalize HumanNet
humannet <- read.table("/scratch/Funcoup5/benchmark/HumanNet-XN.tsv", sep = "\t")
humannet[, 4] <- (humannet[, 3] - min(humannet[, 3])) / (max(humannet[, 3]) - min(humannet[, 3]))
humannet_fil <- humannet[humannet[, 4] >= 0.8, ]

# Load GWAS dataset and filter
gwas <- read.table("/scratch/Funcoup5/benchmark/gwas_catalog_v1.0-associations_e100_r2020-07-14.tsv", header = TRUE, fill = TRUE, sep = "\t")
gwas[, 35] <- as.numeric(substring(gwas[, 1], 1, 2))
new_gwas <- gwas[gwas[, 35] == 20 & as.numeric(gwas[, 28]) < 5e-8, ]
new_gwas[, 35] <- as.numeric(substring(new_gwas[, 4], 1, 4))

# Extract traits and prepare GWAS gene-trait mapping
dat_f <- do.call(rbind, lapply(unique(new_gwas[, 8]), function(trait) {
  ok <- unlist(strsplit(as.vector(new_gwas[new_gwas[, 8] == trait, 18]), ","))
  data.frame(values = trimws(ok), disease = trait)
}))
dat_f <- dat_f[!duplicated(dat_f), ]

# Map Ensembl to Entrez
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
map_ensembl <- getBM(attributes = c('entrezgene_id', 'ensembl_gene_id'), mart = ensembl)
names(dat_f) <- c("ensembl_gene_id", "disease")
translate_diseases <- left_join(dat_f, map_ensembl, by = "ensembl_gene_id") %>% 
  select(entrezgene_id, disease) %>% 
  filter(complete.cases(.)) %>% 
  distinct()

# Load networks
funcoup <- read.delim("/scratch/Funcoup_versions_analysis/fc4.1_human", header = TRUE)[, c(3, 4, 1)]
funcoup <- funcoup[funcoup[, 3] >= 0.8, ]
humannet_genes <- unique(c(humannet_fil[, 1], humannet_fil[, 2]))
funcoup_genes <- unique(c(funcoup[, 1], funcoup[, 2]))
string <- read.delim("/scratch/UNAdrug/string_translated_0.8", header = TRUE)
string_genes <- unique(c(string[, 1], string[, 2]))

# Load gene ID mapping
mapping <- read.delim("/scratch/PathBIX_all/pathbix/data/genetranslations.csv", sep = ",") %>% 
  filter(V4 == "9606")

# Benchmark function
evaluate_network <- function(traits, network_data, gene_column_names, all_network_genes, mapped = TRUE) {
  results <- lapply(traits, function(trait) {
    genes <- translate_diseases$entrezgene_id[translate_diseases$disease == trait]
    original <- genes
    if (mapped) {
      genes <- unique(mapping[mapping[, 1] %in% genes, 2])
    }
    genes <- unique(genes[genes %in% all_network_genes])
    if (length(genes) < 2 || length(original) < 2) return(NULL)
    comb <- combn(genes, 2)
    interactions <- network_data[network_data[[gene_column_names[1]]] %in% genes &
                                   network_data[[gene_column_names[2]]] %in% genes, ]
    data.frame(count = nrow(interactions), total = ncol(comb), trait = trait)
  })
  do.call(rbind, results)
}

# Run benchmarks
traits <- unique(translate_diseases$disease)
funcoup_results <- evaluate_network(traits, funcoup, c(1, 2), funcoup_genes)
humannet_results <- evaluate_network(traits, humannet_fil, c(1, 2), humannet_genes, mapped = FALSE)
string_results <- evaluate_network(traits, string, c(1, 2), string_genes)

# Merge results
colnames(funcoup_results) <- c("Funcoup_links", "total", "disease")
colnames(humannet_results) <- c("Humannet_links", "total", "disease")
colnames(string_results) <- c("String_links", "total", "disease")

merged <- reduce(list(funcoup_results[, c("disease", "Funcoup_links")],
                      string_results[, c("disease", "String_links")],
                      humannet_results[, c("disease", "Humannet_links")]),
                 full_join, by = "disease")
merged[is.na(merged)] <- 0

# Calculate average link recovery
mean_funcoup <- mean(as.numeric(merged$Funcoup_links)) / nrow(funcoup)
mean_string <- mean(as.numeric(merged$String_links)) / nrow(string)
mean_humannet <- mean(as.numeric(merged$Humannet_links)) / nrow(humannet_fil)

# Output results
list(
  funcoup_mean = mean_funcoup,
  string_mean = mean_string,
  humannet_mean = mean_humannet,
  merged_table = merged
)
