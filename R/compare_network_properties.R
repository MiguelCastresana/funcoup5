# Compare link recovery across biological networks.

library(biomaRt)
library(dplyr)
library(purrr)

config <- list(
  gold_dir = Sys.getenv("FUNCOUP_GOLD_DIR", file.path("data", "gold_standards")),
  benchmark_dir = Sys.getenv("FUNCOUP_BENCHMARK_DIR", file.path("data", "benchmark")),
  funcoup_file = Sys.getenv("FUNCOUP_FC4_FILE", file.path("data", "networks", "fc4.1_human")),
  string_file = Sys.getenv("STRING_FILE", file.path("data", "networks", "string_translated_0.8")),
  mapping_file = Sys.getenv("GENE_TRANSLATIONS_FILE", file.path("data", "genetranslations.csv"))
)

required_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing input file: ", path, call. = FALSE)
  }
  path
}

normalize_score <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

gold_files <- list(
  ppi_fc4 = required_file(file.path(config$gold_dir, "ppi_fc4")),
  ppi_fc5 = required_file(file.path(config$gold_dir, "ppi_FC5")),
  complex_fc4 = required_file(file.path(config$gold_dir, "complex_fc4.txt")),
  complex_fc5 = required_file(file.path(config$gold_dir, "complex_FC5"))
)

ppi_fc4 <- read.delim(gold_files$ppi_fc4)
ppi_fc5 <- read.delim(gold_files$ppi_fc5)
common_ids <- intersect(unique(ppi_fc4[, 1]), unique(ppi_fc5[, 1]))

complex_fc4 <- read.delim(gold_files$complex_fc4, header = FALSE, sep = ",")
complex_fc5 <- read.delim(gold_files$complex_fc5, header = FALSE, sep = ",")

humannet <- read.table(required_file(file.path(config$benchmark_dir, "HumanNet-XN.tsv")), sep = "\t")
humannet[, 4] <- normalize_score(humannet[, 3])
humannet_filtered <- humannet[humannet[, 4] >= 0.8, ]

gwas <- read.table(
  required_file(file.path(config$benchmark_dir, "gwas_catalog_v1.0-associations_e100_r2020-07-14.tsv")),
  header = TRUE,
  fill = TRUE,
  sep = "\t"
)
gwas[, 35] <- as.numeric(substring(gwas[, 1], 1, 2))
gwas_filtered <- gwas[gwas[, 35] == 20 & as.numeric(gwas[, 28]) < 5e-8, ]
gwas_filtered[, 35] <- as.numeric(substring(gwas_filtered[, 4], 1, 4))

gwas_genes <- do.call(rbind, lapply(unique(gwas_filtered[, 8]), function(trait) {
  genes <- unlist(strsplit(as.vector(gwas_filtered[gwas_filtered[, 8] == trait, 18]), ","))
  data.frame(values = trimws(genes), disease = trait)
}))
gwas_genes <- gwas_genes[!duplicated(gwas_genes), ]

ensembl <- biomaRt::useMart("ensembl")
ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)
ensembl_map <- biomaRt::getBM(attributes = c("entrezgene_id", "ensembl_gene_id"), mart = ensembl)

names(gwas_genes) <- c("ensembl_gene_id", "disease")
translated_diseases <- gwas_genes %>%
  left_join(ensembl_map, by = "ensembl_gene_id") %>%
  select(entrezgene_id, disease) %>%
  filter(complete.cases(.)) %>%
  distinct()

funcoup <- read.delim(required_file(config$funcoup_file), header = TRUE)[, c(3, 4, 1)]
funcoup <- funcoup[funcoup[, 3] >= 0.8, ]

string <- read.delim(required_file(config$string_file), header = TRUE)
mapping <- read.delim(required_file(config$mapping_file), sep = ",") %>%
  filter(V4 == "9606")

network_genes <- list(
  funcoup = unique(c(funcoup[, 1], funcoup[, 2])),
  string = unique(c(string[, 1], string[, 2])),
  humannet = unique(c(humannet_filtered[, 1], humannet_filtered[, 2]))
)

evaluate_network <- function(traits, network_data, gene_columns, all_network_genes, mapped = TRUE) {
  results <- lapply(traits, function(trait) {
    genes <- translated_diseases$entrezgene_id[translated_diseases$disease == trait]
    original_genes <- genes

    if (mapped) {
      genes <- unique(mapping[mapping[, 1] %in% genes, 2])
    }

    genes <- unique(genes[genes %in% all_network_genes])
    if (length(genes) < 2 || length(original_genes) < 2) {
      return(NULL)
    }

    possible_pairs <- utils::combn(genes, 2)
    interactions <- network_data[
      network_data[[gene_columns[1]]] %in% genes &
        network_data[[gene_columns[2]]] %in% genes,
    ]

    data.frame(count = nrow(interactions), total = ncol(possible_pairs), trait = trait)
  })

  do.call(rbind, results)
}

traits <- unique(translated_diseases$disease)
funcoup_results <- evaluate_network(traits, funcoup, c(1, 2), network_genes$funcoup)
humannet_results <- evaluate_network(traits, humannet_filtered, c(1, 2), network_genes$humannet, mapped = FALSE)
string_results <- evaluate_network(traits, string, c(1, 2), network_genes$string)

colnames(funcoup_results) <- c("FunCoup_links", "total", "disease")
colnames(humannet_results) <- c("HumanNet_links", "total", "disease")
colnames(string_results) <- c("STRING_links", "total", "disease")

merged <- purrr::reduce(
  list(
    funcoup_results[, c("disease", "FunCoup_links")],
    string_results[, c("disease", "STRING_links")],
    humannet_results[, c("disease", "HumanNet_links")]
  ),
  full_join,
  by = "disease"
)
merged[is.na(merged)] <- 0

list(
  common_gold_standard_ids = common_ids,
  complex_sets = list(fc4 = complex_fc4, fc5 = complex_fc5),
  funcoup_mean = mean(as.numeric(merged$FunCoup_links)) / nrow(funcoup),
  string_mean = mean(as.numeric(merged$STRING_links)) / nrow(string),
  humannet_mean = mean(as.numeric(merged$HumanNet_links)) / nrow(humannet_filtered),
  merged_table = merged
)
