# Random walk with restart across FunCoup, STRING, and HumanNet-style networks.

library(dnet)
library(PRROC)
library(igraph)

load_network <- function(path, sep = "\t", header = TRUE) {
  network <- read.delim(path, sep = sep, header = header, stringsAsFactors = FALSE)
  network <- network[, 1:2, drop = FALSE]
  network <- stats::na.omit(network)
  unique(network)
}

run_affinity <- function(graph, seeds, end_points, genes = igraph::V(graph)$name) {
  if (is.null(genes) || anyNA(genes)) {
    stop("Graph vertices must have names or `genes` must be provided.", call. = FALSE)
  }

  seeds <- intersect(seeds, genes)
  end_points <- intersect(end_points, genes)
  if (length(seeds) == 0 || length(end_points) == 0) {
    return(list(pr_auc = NA_real_, roc_auc = NA_real_))
  }

  seed_vector <- as.integer(genes %in% seeds)
  seed_data <- data.frame(seed = seed_vector, row.names = genes)

  affinity <- dnet::dRWR(
    g = graph,
    normalise = "laplacian",
    setSeeds = seed_data,
    restart = 0.75,
    parallel = TRUE
  )

  scores <- data.frame(score = affinity[, 1], gene = genes)
  scores <- scores[order(-scores$score), ]

  true_scores <- scores$score[scores$gene %in% end_points]
  false_scores <- scores$score[!(scores$gene %in% end_points)]

  if (length(true_scores) == 0 || length(false_scores) == 0) {
    return(list(pr_auc = NA_real_, roc_auc = NA_real_))
  }

  pr <- PRROC::pr.curve(scores.class0 = true_scores, scores.class1 = false_scores)
  roc <- PRROC::roc.curve(scores.class0 = true_scores, scores.class1 = false_scores)

  list(pr_auc = pr$auc.integral, roc_auc = roc$auc)
}

create_splits <- function(kegg, n = 30, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  groups <- unique(as.vector(kegg[, 2]))
  groups <- groups[vapply(groups, function(group) {
    length(unique(as.vector(kegg[kegg[, 2] == group, 1]))) >= 2
  }, logical(1))]

  if (length(groups) == 0) {
    stop("No KEGG groups contain at least two genes.", call. = FALSE)
  }

  data_a <- vector("list", n)
  data_b <- vector("list", n)

  for (split in seq_len(n)) {
    split_a <- data.frame()
    split_b <- data.frame()

    for (group in groups) {
      genes <- unique(as.vector(kegg[kegg[, 2] == group, 1]))
      genes_a <- sample(genes, floor(length(genes) / 2))
      genes_b <- setdiff(genes, genes_a)

      split_a <- rbind(split_a, data.frame(ensembl_gene_id = genes_a, disease = group))
      split_b <- rbind(split_b, data.frame(ensembl_gene_id = genes_b, disease = group))
    }

    data_a[[split]] <- split_a
    data_b[[split]] <- split_b
  }

  list(A = data_a, B = data_b)
}

benchmark_network <- function(graph_data, kegg, splits, output_prefix, n_networks = 30) {
  if (nrow(graph_data) == 0) {
    stop("`graph_data` must contain at least one edge.", call. = FALSE)
  }
  if (missing(splits) || !all(c("A", "B") %in% names(splits))) {
    stop("`splits` must be a list with A and B entries.", call. = FALSE)
  }

  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  results_pr <- vector("list", n_networks)
  results_roc <- vector("list", n_networks)

  for (network_id in seq_len(n_networks)) {
    graph <- igraph::graph_from_data_frame(graph_data[, 1:2], directed = FALSE)
    graph <- igraph::rewire(graph, with = igraph::keeping_degseq(niter = nrow(graph_data) * 100))
    genes <- igraph::V(graph)$name

    pr_values <- numeric()
    roc_values <- numeric()

    for (i in seq_along(splits$A)) {
      part_a <- split(splits$A[[i]], splits$A[[i]][, 2])
      part_b <- split(splits$B[[i]], splits$B[[i]][, 2])

      for (group in intersect(names(part_a), names(part_b))) {
        result <- run_affinity(graph, part_a[[group]][, 1], part_b[[group]][, 1], genes)
        pr_values <- c(pr_values, result$pr_auc)
        roc_values <- c(roc_values, result$roc_auc)
      }
    }

    results_pr[[network_id]] <- pr_values
    results_roc[[network_id]] <- roc_values

    saveRDS(pr_values, file = paste0(output_prefix, "_pr_network_", network_id, ".rds"))
    saveRDS(roc_values, file = paste0(output_prefix, "_roc_network_", network_id, ".rds"))
  }

  invisible(list(pr = results_pr, roc = results_roc))
}
