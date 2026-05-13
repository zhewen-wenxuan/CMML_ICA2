############################################################
## CMML3 ICA2 - Miniproject 3
## scDesign3 simulation + trajectory inference benchmarking
##
## Project logic:
## 1. Load scDesign3 tutorial pancreatic endocrinogenesis data
## 2. Simulate baseline synthetic scRNA-seq data using scDesign3
## 3. Validate reference vs synthetic similarity
## 4. Create three benchmark scenarios:
##    S1 baseline, S2 low-depth/high-noise, S3 imbalanced trajectory
## 5. Benchmark Slingshot, Monocle3 and TSCAN
## 6. Generate main and supplementary figures
############################################################


############################################################
## 0. Settings
############################################################

set.seed(123)

project_dir <- "CMML3_ICA2_scDesign3_trajectory"
dir.create(project_dir, showWarnings = FALSE)
dir.create(file.path(project_dir, "data"), showWarnings = FALSE)
dir.create(file.path(project_dir, "data", "processed"), showWarnings = FALSE)
dir.create(file.path(project_dir, "results"), showWarnings = FALSE)
dir.create(file.path(project_dir, "results", "objects"), showWarnings = FALSE)
dir.create(file.path(project_dir, "results", "metrics"), showWarnings = FALSE)
dir.create(file.path(project_dir, "results", "figures"), showWarnings = FALSE)
dir.create(file.path(project_dir, "results", "supplementary"), showWarnings = FALSE)

## Keep this small first. The course quickstart also uses 100 genes.
n_genes_use <- 100

## Number of PCs used for trajectory methods
n_pcs_use <- 20

## Number of k-means clusters used for Slingshot input
k_clusters_use <- 6

## Low-depth scenario: keep 25% of counts (harder)
downsample_prob <- 0.25

## Imbalanced scenario: keep only 15% of middle pseudotime cells (harder)
middle_keep_prob <- 0.15

## Install missing packages automatically?
INSTALL_MISSING_PACKAGES <- TRUE


############################################################
## 1. Package installation and loading
############################################################

install_if_missing_cran <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_missing_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      BiocManager::install(p, update = FALSE, ask = FALSE)
    }
  }
}

if (INSTALL_MISSING_PACKAGES) {
  install_if_missing_cran(c(
    "ggplot2", "dplyr", "tidyr", "readr", "patchwork",
    "mclust", "uwot", "Matrix", "remotes"
  ))
  
  install_if_missing_bioc(c(
    "SingleCellExperiment", "SummarizedExperiment",
    "S4Vectors", "scDesign3", "slingshot", "TSCAN"
  ))
  
  ## Monocle3 installation is sometimes the most fragile part.
  ## Try Bioconductor first; if not available, try GitHub.
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    tryCatch({
      BiocManager::install("monocle3", update = FALSE, ask = FALSE)
    }, error = function(e) {
      message("Bioconductor monocle3 installation failed. Trying GitHub...")
    })
  }
  
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    remotes::install_github("cole-trapnell-lab/monocle3")
  }
}


install.packages("BiocManager")
BiocManager::install(version = "3.20", ask = FALSE)
## 关键：让 remotes 能找到 Bioconductor 包
options(repos = BiocManager::repositories())
## 先单独安装 monocle3 缺的依赖
BiocManager::install("batchelor", ask = FALSE, update = FALSE)
## 再装 monocle3
install.packages("remotes")
remotes::install_github(
  "cole-trapnell-lab/monocle3",
  dependencies = TRUE,
  upgrade = "never"
)

library(monocle3)



library(SingleCellExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(scDesign3)
library(slingshot)
library(TSCAN)
library(monocle3)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(mclust)
library(uwot)
library(Matrix)

theme_set(theme_bw(base_size = 11))


############################################################
## 2. Helper functions
############################################################

scale01 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

safe_num <- function(x) {
  as.numeric(as.character(x))
}

make_logcounts <- function(sce) {
  logcounts(sce) <- log1p(as.matrix(counts(sce)))
  sce
}

calc_cell_qc <- function(sce, dataset_name) {
  mat <- as.matrix(counts(sce))
  data.frame(
    cell = colnames(mat),
    dataset = dataset_name,
    library_size = colSums(mat),
    detected_genes = colSums(mat > 0),
    zero_proportion = colMeans(mat == 0),
    pseudotime = scale01(safe_num(colData(sce)$pseudotime)),
    stringsAsFactors = FALSE
  )
}

calc_gene_mean_var <- function(sce, dataset_name) {
  mat <- as.matrix(counts(sce))
  data.frame(
    gene = rownames(mat),
    dataset = dataset_name,
    mean_count = rowMeans(mat),
    variance_count = apply(mat, 1, var),
    stringsAsFactors = FALSE
  )
}

create_sce_from_counts <- function(count_mat, coldata_df, scenario_name) {
  count_mat <- round(as.matrix(count_mat))
  count_mat[count_mat < 0] <- 0
  
  if (is.null(rownames(count_mat))) {
    rownames(count_mat) <- paste0("Gene", seq_len(nrow(count_mat)))
  }
  if (is.null(colnames(count_mat))) {
    colnames(count_mat) <- paste0("Cell", seq_len(ncol(count_mat)))
  }
  
  coldata_df <- as.data.frame(coldata_df)
  rownames(coldata_df) <- colnames(count_mat)
  coldata_df$scenario <- scenario_name
  
  if (!"pseudotime" %in% colnames(coldata_df)) {
    stop("pseudotime is missing from colData. This project needs true pseudotime.")
  }
  coldata_df$pseudotime <- scale01(safe_num(coldata_df$pseudotime))
  
  sce <- SingleCellExperiment(
    assays = list(counts = count_mat),
    colData = S4Vectors::DataFrame(coldata_df)
  )
  logcounts(sce) <- log1p(counts(sce))
  sce
}

make_pca_umap_clusters <- function(sce, n_pcs = 20, k_clusters = 6) {
  sce <- make_logcounts(sce)
  
  mat <- t(as.matrix(logcounts(sce)))   # cells x genes
  vars <- apply(mat, 2, var, na.rm = TRUE)
  mat <- mat[, vars > 0, drop = FALSE]
  
  if (ncol(mat) < 2) {
    stop("Too few variable genes for PCA.")
  }
  
  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  n_pcs_real <- min(n_pcs, ncol(pca$x))
  pcs <- pca$x[, seq_len(n_pcs_real), drop = FALSE]
  
  reducedDim(sce, "PCA") <- pcs
  
  set.seed(123)
  um <- uwot::umap(
    pcs,
    n_neighbors = min(30, nrow(pcs) - 1),
    min_dist = 0.3,
    metric = "cosine"
  )
  rownames(um) <- colnames(sce)
  colnames(um) <- c("UMAP1", "UMAP2")
  reducedDim(sce, "UMAP") <- um
  
  set.seed(123)
  k_real <- min(k_clusters, ncol(sce) - 1)
  km <- kmeans(pcs[, seq_len(min(10, ncol(pcs))), drop = FALSE],
               centers = k_real,
               nstart = 20)
  colData(sce)$cluster_kmeans <- factor(km$cluster)
  
  sce
}

bin_pseudotime_three <- function(x) {
  x <- scale01(x)
  q <- quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  
  ## In rare cases duplicated breaks can happen; use rank fallback.
  if (length(unique(q)) < 4) {
    x_rank <- rank(x, ties.method = "average")
    q <- quantile(x_rank, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
    return(cut(x_rank, breaks = q, include.lowest = TRUE,
               labels = c("early", "middle", "late")))
  }
  
  cut(x, breaks = q, include.lowest = TRUE,
      labels = c("early", "middle", "late"))
}

orient_pseudotime <- function(true_pt, pred_pt) {
  true_pt <- scale01(true_pt)
  pred_pt <- scale01(pred_pt)
  
  ok <- is.finite(true_pt) & is.finite(pred_pt)
  if (sum(ok) < 5) return(pred_pt)
  
  rho <- suppressWarnings(cor(true_pt[ok], pred_pt[ok], method = "spearman"))
  if (is.finite(rho) && rho < 0) {
    pred_pt <- 1 - pred_pt
  }
  pred_pt
}

calc_metrics_one <- function(true_pt, pred_pt) {
  true_pt <- scale01(true_pt)
  pred_pt <- orient_pseudotime(true_pt, pred_pt)
  
  ok <- is.finite(true_pt) & is.finite(pred_pt)
  true_ok <- true_pt[ok]
  pred_ok <- pred_pt[ok]
  
  if (length(true_ok) < 5) {
    return(data.frame(
      spearman = NA_real_,
      kendall = NA_real_,
      rmse = NA_real_,
      ari_three_state = NA_real_
    ))
  }
  
  spearman <- suppressWarnings(cor(true_ok, pred_ok, method = "spearman"))
  kendall <- suppressWarnings(cor(true_ok, pred_ok, method = "kendall"))
  rmse <- sqrt(mean((true_ok - pred_ok)^2))
  
  true_bin <- bin_pseudotime_three(true_ok)
  pred_bin <- bin_pseudotime_three(pred_ok)
  ari <- mclust::adjustedRandIndex(true_bin, pred_bin)
  
  data.frame(
    spearman = spearman,
    kendall = kendall,
    rmse = rmse,
    ari_three_state = ari
  )
}


############################################################
## 3. Load scDesign3 tutorial data
############################################################

message("Loading scDesign3 tutorial dataset...")

example_sce <- readRDS("/Users/xuanzhewen/Downloads/scDesign3_pancreas_example_sce.rds")

message("Original dataset:")
print(example_sce)

## Use first 100 genes, following the course quickstart
example_sce <- example_sce[seq_len(min(n_genes_use, nrow(example_sce))), ]

## Make sure a counts assay exists
if (!"counts" %in% assayNames(example_sce)) {
  if ("X" %in% assayNames(example_sce)) {
    assay(example_sce, "counts") <- round(as.matrix(assay(example_sce, "X")))
  } else {
    stop("No counts assay found. Please inspect assayNames(example_sce).")
  }
}

counts(example_sce) <- round(as.matrix(counts(example_sce)))

if (!"pseudotime" %in% colnames(colData(example_sce))) {
  stop("The reference dataset does not contain colData$pseudotime.")
}

colData(example_sce)$pseudotime <- scale01(safe_num(colData(example_sce)$pseudotime))
example_sce <- make_logcounts(example_sce)

saveRDS(example_sce, file.path(project_dir, "data", "processed", "reference_example_sce.rds"))


############################################################
## 4. Run scDesign3 baseline simulation
############################################################

message("Running scDesign3 simulation... This may take several minutes.")

## The course slides use celltype = "cell_type", but the object may not always
## contain this exact column. Use it only if present.
celltype_arg <- if ("cell_type" %in% colnames(colData(example_sce))) {
  "cell_type"
} else {
  NULL
}

set.seed(123)
example_simu <- scdesign3(
  sce = example_sce,
  assay_use = "counts",
  celltype = celltype_arg,
  pseudotime = "pseudotime",
  spatial = NULL,
  other_covariates = NULL,
  mu_formula = "s(pseudotime, k = 10, bs = 'cr')",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = FALSE
)

sim_counts <- example_simu$new_count
sim_counts <- round(as.matrix(sim_counts))

## If scDesign3 did not generate new covariates, use the original covariates.
## This happens when ncell is not changed.
if (is.null(example_simu$new_covariate)) {
  sim_coldata <- as.data.frame(colData(example_sce))
} else {
  sim_coldata <- as.data.frame(example_simu$new_covariate)
}

rownames(sim_coldata) <- colnames(sim_counts)

sce_s1 <- create_sce_from_counts(
  count_mat = sim_counts,
  coldata_df = sim_coldata,
  scenario_name = "S1_baseline"
)

sce_s1 <- make_pca_umap_clusters(sce_s1, n_pcs = n_pcs_use, k_clusters = k_clusters_use)

saveRDS(example_simu, file.path(project_dir, "results", "objects", "scdesign3_output_list.rds"))
saveRDS(sce_s1, file.path(project_dir, "results", "objects", "sce_S1_baseline.rds"))


############################################################
## 5. Create S2 low-depth / high-noise scenario
############################################################

message("Creating S2 low-depth/high-noise scenario...")

set.seed(123)
s1_counts <- as.matrix(counts(sce_s1))

s2_counts <- matrix(
  rbinom(
    n = length(s1_counts),
    size = as.vector(s1_counts),
    prob = downsample_prob
  ),
  nrow = nrow(s1_counts),
  ncol = ncol(s1_counts),
  byrow = FALSE
)
rownames(s2_counts) <- rownames(s1_counts)
colnames(s2_counts) <- colnames(s1_counts)

s2_coldata <- as.data.frame(colData(sce_s1))
sce_s2 <- create_sce_from_counts(
  count_mat = s2_counts,
  coldata_df = s2_coldata,
  scenario_name = "S2_low_depth"
)
sce_s2 <- make_pca_umap_clusters(sce_s2, n_pcs = n_pcs_use, k_clusters = k_clusters_use)

saveRDS(sce_s2, file.path(project_dir, "results", "objects", "sce_S2_low_depth.rds"))


############################################################
## 6. Create S3 imbalanced trajectory scenario
############################################################

message("Creating S3 imbalanced trajectory scenario...")

set.seed(123)
pt_s1 <- scale01(safe_num(colData(sce_s1)$pseudotime))
pt_bin <- bin_pseudotime_three(pt_s1)

keep_cells <- rep(TRUE, length(pt_bin))
keep_cells[pt_bin == "middle"] <- runif(sum(pt_bin == "middle")) < middle_keep_prob
keep_cells <- keep_cells & !is.na(pt_bin)

sce_s3 <- sce_s1[, keep_cells]
colData(sce_s3)$scenario <- "S3_imbalanced_middle"
sce_s3 <- make_pca_umap_clusters(sce_s3, n_pcs = n_pcs_use, k_clusters = k_clusters_use)

saveRDS(sce_s3, file.path(project_dir, "results", "objects", "sce_S3_imbalanced_middle.rds"))


############################################################
## 7. Reference and scenario QC
############################################################

message("Calculating QC summaries...")

example_sce_for_qc <- make_pca_umap_clusters(
  example_sce,
  n_pcs = n_pcs_use,
  k_clusters = k_clusters_use
)

qc_reference <- calc_cell_qc(example_sce_for_qc, "Reference")
qc_s1 <- calc_cell_qc(sce_s1, "S1_baseline")
qc_s2 <- calc_cell_qc(sce_s2, "S2_low_depth")
qc_s3 <- calc_cell_qc(sce_s3, "S3_imbalanced_middle")

qc_all <- bind_rows(qc_reference, qc_s1, qc_s2, qc_s3)

write_csv(qc_all, file.path(project_dir, "results", "metrics", "cell_qc_all.csv"))

gene_mv_reference <- calc_gene_mean_var(example_sce_for_qc, "Reference")
gene_mv_s1 <- calc_gene_mean_var(sce_s1, "S1_baseline")
gene_mv_all <- bind_rows(gene_mv_reference, gene_mv_s1)

write_csv(gene_mv_all, file.path(project_dir, "results", "metrics", "gene_mean_variance_reference_vs_s1.csv"))


############################################################
## 8. Trajectory inference wrappers
############################################################

run_slingshot_pt <- function(sce) {
  sce <- make_pca_umap_clusters(sce, n_pcs = n_pcs_use, k_clusters = k_clusters_use)
  
  true_pt <- scale01(safe_num(colData(sce)$pseudotime))
  cl <- colData(sce)$cluster_kmeans
  
  mean_pt_by_cluster <- tapply(true_pt, cl, mean, na.rm = TRUE)
  start_cluster <- names(which.min(mean_pt_by_cluster))
  
  sce_sl <- slingshot(
    sce,
    clusterLabels = "cluster_kmeans",
    reducedDim = "PCA",
    start.clus = start_cluster
  )
  
  pt_mat <- slingPseudotime(sce_sl)
  
  ## If multiple lineages are returned, choose the lineage with most finite cells.
  finite_counts <- colSums(is.finite(pt_mat))
  lineage_use <- which.max(finite_counts)
  
  pred <- pt_mat[, lineage_use]
  names(pred) <- colnames(sce)
  pred
}

run_monocle3_pt <- function(sce) {
  count_mat <- as.matrix(counts(sce))
  cell_meta <- as.data.frame(colData(sce))
  gene_meta <- data.frame(
    gene_short_name = rownames(count_mat),
    row.names = rownames(count_mat)
  )
  
  cds <- new_cell_data_set(
    expression_data = count_mat,
    cell_metadata = cell_meta,
    gene_metadata = gene_meta
  )
  
  cds <- preprocess_cds(cds, num_dim = min(n_pcs_use, nrow(count_mat) - 1))
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = FALSE)
  
  true_pt <- scale01(safe_num(colData(sce)$pseudotime))
  root_cells <- colnames(sce)[true_pt <= quantile(true_pt, 0.05, na.rm = TRUE)]
  
  cds <- order_cells(cds, root_cells = root_cells)
  pred <- pseudotime(cds)
  pred <- as.numeric(pred[colnames(sce)])
  names(pred) <- colnames(sce)
  pred
}

run_tscan_pt <- function(sce) {
  expr <- as.matrix(logcounts(sce))  # genes x cells
  
  ## TSCAN::preprocess takes genes x cells expression matrix.
  ## takelog = FALSE because logcounts are already log1p transformed.
  ## 加上 cvcutoff = 0 和 minexpr_percent = 0 停止过滤基因，防止基因太少(100)被全部删掉导致报错
  proc <- TSCAN::preprocess(expr, takelog = FALSE, cvcutoff = 0, minexpr_percent = 0)
  mcl <- TSCAN::exprmclust(proc)
  ord <- TSCAN::TSCANorder(mcl)
  
  pred <- rep(NA_real_, ncol(sce))
  names(pred) <- colnames(sce)
  
  ord <- intersect(ord, colnames(sce))
  pred[ord] <- seq_along(ord)
  
  pred
}

run_all_methods_one_scenario <- function(sce, scenario_name) {
  message("Running trajectory inference for: ", scenario_name)
  
  true_pt <- scale01(safe_num(colData(sce)$pseudotime))
  names(true_pt) <- colnames(sce)
  
  pred_list <- list()
  
  pred_list$Slingshot <- tryCatch({
    run_slingshot_pt(sce)
  }, error = function(e) {
    message("Slingshot failed for ", scenario_name, ": ", e$message)
    x <- rep(NA_real_, ncol(sce)); names(x) <- colnames(sce); x
  })
  
  pred_list$Monocle3 <- tryCatch({
    run_monocle3_pt(sce)
  }, error = function(e) {
    message("Monocle3 failed for ", scenario_name, ": ", e$message)
    x <- rep(NA_real_, ncol(sce)); names(x) <- colnames(sce); x
  })
  
  pred_list$TSCAN <- tryCatch({
    run_tscan_pt(sce)
  }, error = function(e) {
    message("TSCAN failed for ", scenario_name, ": ", e$message)
    x <- rep(NA_real_, ncol(sce)); names(x) <- colnames(sce); x
  })
  
  pred_df <- bind_rows(lapply(names(pred_list), function(method_name) {
    pred <- pred_list[[method_name]]
    pred <- pred[colnames(sce)]
    pred_aligned <- orient_pseudotime(true_pt, pred)
    
    data.frame(
      cell = colnames(sce),
      scenario = scenario_name,
      method = method_name,
      true_pseudotime = as.numeric(true_pt),
      inferred_pseudotime_raw = as.numeric(scale01(pred)),
      inferred_pseudotime = as.numeric(pred_aligned),
      stringsAsFactors = FALSE
    )
  }))
  
  metric_df <- bind_rows(lapply(names(pred_list), function(method_name) {
    pred <- pred_list[[method_name]]
    pred <- pred[colnames(sce)]
    m <- calc_metrics_one(true_pt, pred)
    data.frame(
      scenario = scenario_name,
      method = method_name,
      m,
      stringsAsFactors = FALSE
    )
  }))
  
  list(pred = pred_df, metrics = metric_df)
}


############################################################
## 9. Run Slingshot / Monocle3 / TSCAN on three scenarios
############################################################

scenario_list <- list(
  S1_baseline = sce_s1,
  S2_low_depth = sce_s2,
  S3_imbalanced_middle = sce_s3
)

ti_results <- lapply(names(scenario_list), function(scn) {
  run_all_methods_one_scenario(scenario_list[[scn]], scn)
})
names(ti_results) <- names(scenario_list)

pred_all <- bind_rows(lapply(ti_results, function(x) x$pred))
metrics_all <- bind_rows(lapply(ti_results, function(x) x$metrics))

write_csv(pred_all, file.path(project_dir, "results", "metrics", "pseudotime_predictions_all.csv"))
write_csv(metrics_all, file.path(project_dir, "results", "metrics", "trajectory_metrics_all.csv"))

saveRDS(ti_results, file.path(project_dir, "results", "objects", "trajectory_inference_results.rds"))


############################################################
## 10. Main Figure 1: simulation design and validation
############################################################

message("Generating Figure 1...")

## Figure 1A: workflow schematic
workflow_df <- data.frame(
  x = c(1, 2, 3, 4, 5),
  y = rep(1, 5),
  label = c(
    "Reference\nscRNA-seq\nwith pseudotime",
    "scDesign3\nbaseline\nsimulation",
    "Three\nsimulation\nscenarios",
    "Slingshot\nMonocle3\nTSCAN",
    "Pseudotime\nrecovery\nmetrics"
  )
)

p1a <- ggplot(workflow_df, aes(x, y)) +
  geom_label(aes(label = label), size = 3.2, label.size = 0.3) +
  geom_segment(
    data = data.frame(x = 1:4, xend = 2:5, y = 1, yend = 1),
    aes(x = x + 0.35, xend = xend - 0.35, y = y, yend = yend),
    arrow = arrow(length = unit(0.18, "cm")),
    inherit.aes = FALSE
  ) +
  xlim(0.5, 5.5) +
  ylim(0.7, 1.3) +
  labs(title = "A  Simulation-based benchmarking workflow") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

## Figure 1B: reference vs S1 UMAP
umap_ref <- as.data.frame(reducedDim(example_sce_for_qc, "UMAP"))
colnames(umap_ref) <- c("UMAP1", "UMAP2")
umap_ref$dataset <- "Reference"
umap_ref$pseudotime <- scale01(safe_num(colData(example_sce_for_qc)$pseudotime))

umap_s1 <- as.data.frame(reducedDim(sce_s1, "UMAP"))
colnames(umap_s1) <- c("UMAP1", "UMAP2")
umap_s1$dataset <- "S1 baseline synthetic"
umap_s1$pseudotime <- scale01(safe_num(colData(sce_s1)$pseudotime))

umap_ref_s1 <- bind_rows(umap_ref, umap_s1)

p1b <- ggplot(umap_ref_s1, aes(UMAP1, UMAP2, colour = pseudotime)) +
  geom_point(size = 0.35, alpha = 0.8) +
  facet_wrap(~ dataset, nrow = 1) +
  scale_colour_viridis_c(option = "viridis") +
  labs(
    title = "B  Reference and baseline synthetic data",
    colour = "True\npseudotime"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90")
  )

## Figure 1C: QC comparison reference vs S1
qc_ref_s1_long <- qc_all %>%
  filter(dataset %in% c("Reference", "S1_baseline")) %>%
  mutate(
    log10_library_size = log10(library_size + 1)
  ) %>%
  select(dataset, log10_library_size, detected_genes, zero_proportion) %>%
  pivot_longer(
    cols = c(log10_library_size, detected_genes, zero_proportion),
    names_to = "metric",
    values_to = "value"
  )

p1c <- ggplot(qc_ref_s1_long, aes(dataset, value)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, linewidth = 0.2) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(
    title = "C  Cell-level QC similarity",
    x = NULL,
    y = "Value"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "grey90")
  )

## Figure 1D: scenario design summary
qc_scenarios <- qc_all %>%
  filter(dataset %in% c("S1_baseline", "S2_low_depth", "S3_imbalanced_middle")) %>%
  mutate(pt_bin = bin_pseudotime_three(pseudotime))

scenario_prop <- qc_scenarios %>%
  dplyr::count(dataset, pt_bin) %>%
  dplyr::group_by(dataset) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

p1d <- ggplot(scenario_prop, aes(dataset, prop, fill = pt_bin)) +
  geom_col(position = "stack", colour = "white", linewidth = 0.2) +
  labs(
    title = "D  Pseudotime-state proportions across scenarios",
    x = NULL,
    y = "Proportion",
    fill = "True\nstate"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

fig1 <- p1a / p1b / p1c / p1d +
  plot_layout(heights = c(0.7, 1.4, 1.1, 1.1))

ggsave(
  filename = file.path(project_dir, "results", "figures", "Figure1_simulation_validation.pdf"),
  plot = fig1,
  width = 11,
  height = 12
)

ggsave(
  filename = file.path(project_dir, "results", "figures", "Figure1_simulation_validation.png"),
  plot = fig1,
  width = 11,
  height = 12,
  dpi = 300
)


############################################################
## 11. Main Figure 2: trajectory benchmarking
############################################################

message("Generating Figure 2...")

metrics_long <- metrics_all %>%
  pivot_longer(
    cols = c(spearman, kendall, rmse, ari_three_state),
    names_to = "metric",
    values_to = "value"
  )

## Figure 2A: Spearman heatmap
p2a_df <- metrics_all %>%
  select(scenario, method, spearman)

p2a <- ggplot(p2a_df, aes(scenario, method, fill = spearman)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", spearman)), size = 3.2) +
  scale_fill_viridis_c(option = "viridis", limits = c(0, 1), na.value = "grey80") +
  labs(
    title = "A  Spearman correlation with true pseudotime",
    x = NULL,
    y = NULL,
    fill = "Spearman"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

## Figure 2B: RMSE bar plot
p2b <- ggplot(metrics_all, aes(scenario, rmse, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  labs(
    title = "B  Scaled pseudotime RMSE",
    x = NULL,
    y = "RMSE",
    fill = "Method"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

## Figure 2C: representative scatter plot for hardest scenario
hard_scenario <- "S3_imbalanced_middle"

p2c_df <- pred_all %>%
  filter(scenario == hard_scenario) %>%
  filter(is.finite(inferred_pseudotime))

p2c <- ggplot(p2c_df, aes(true_pseudotime, inferred_pseudotime)) +
  geom_point(size = 0.35, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ method, nrow = 1) +
  labs(
    title = "C  True vs inferred pseudotime in imbalanced scenario",
    x = "True pseudotime",
    y = "Inferred pseudotime"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90")
  )

## Figure 2D: ARI for early/middle/late recovery
p2d <- ggplot(metrics_all, aes(scenario, ari_three_state, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  labs(
    title = "D  Early/middle/late state recovery",
    x = NULL,
    y = "ARI",
    fill = "Method"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

fig2 <- (p2a | p2b) / (p2c / p2d) +
  plot_layout(heights = c(1, 1.4))

ggsave(
  filename = file.path(project_dir, "results", "figures", "Figure2_trajectory_benchmarking.pdf"),
  plot = fig2,
  width = 12,
  height = 10
)

ggsave(
  filename = file.path(project_dir, "results", "figures", "Figure2_trajectory_benchmarking.png"),
  plot = fig2,
  width = 12,
  height = 10,
  dpi = 300
)


############################################################
## 12. Supplementary Figure 1: QC across all scenarios
############################################################

message("Generating supplementary figures...")

qc_all_long <- qc_all %>%
  mutate(log10_library_size = log10(library_size + 1)) %>%
  select(dataset, log10_library_size, detected_genes, zero_proportion) %>%
  pivot_longer(
    cols = c(log10_library_size, detected_genes, zero_proportion),
    names_to = "metric",
    values_to = "value"
  )

sup1 <- ggplot(qc_all_long, aes(dataset, value)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, linewidth = 0.2) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(
    title = "Supplementary Figure 1. Cell-level QC across reference and simulated scenarios",
    x = NULL,
    y = "Value"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  filename = file.path(project_dir, "results", "supplementary", "Supplementary_Figure1_QC_all_scenarios.pdf"),
  plot = sup1,
  width = 11,
  height = 4
)

ggsave(
  filename = file.path(project_dir, "results", "supplementary", "Supplementary_Figure1_QC_all_scenarios.png"),
  plot = sup1,
  width = 11,
  height = 4,
  dpi = 300
)


############################################################
## 13. Supplementary Figure 2: all true vs inferred pseudotime
############################################################

sup2 <- pred_all %>%
  filter(is.finite(inferred_pseudotime)) %>%
  ggplot(aes(true_pseudotime, inferred_pseudotime)) +
  geom_point(size = 0.25, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_grid(method ~ scenario) +
  labs(
    title = "Supplementary Figure 2. True vs inferred pseudotime for all methods and scenarios",
    x = "True pseudotime",
    y = "Inferred pseudotime"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90")
  )

ggsave(
  filename = file.path(project_dir, "results", "supplementary", "Supplementary_Figure2_all_pseudotime_scatter.pdf"),
  plot = sup2,
  width = 12,
  height = 8
)

ggsave(
  filename = file.path(project_dir, "results", "supplementary", "Supplementary_Figure2_all_pseudotime_scatter.png"),
  plot = sup2,
  width = 12,
  height = 8,
  dpi = 300
)


############################################################
## 14. Supplementary Figure 3: UMAPs of three scenarios
############################################################

get_umap_scenario_df <- function(sce, scenario_name) {
  sce <- make_pca_umap_clusters(sce, n_pcs = n_pcs_use, k_clusters = k_clusters_use)
  um <- as.data.frame(reducedDim(sce, "UMAP"))
  colnames(um) <- c("UMAP1", "UMAP2")
  um$scenario <- scenario_name
  um$pseudotime <- scale01(safe_num(colData(sce)$pseudotime))
  um
}

umap_scenarios <- bind_rows(
  get_umap_scenario_df(sce_s1, "S1 baseline"),
  get_umap_scenario_df(sce_s2, "S2 low-depth"),
  get_umap_scenario_df(sce_s3, "S3 imbalanced")
)

sup3 <- ggplot(umap_scenarios, aes(UMAP1, UMAP2, colour = pseudotime)) +
  geom_point(size = 0.35, alpha = 0.8) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_colour_viridis_c(option = "viridis") +
  labs(
    title = "Supplementary Figure 3. UMAP visualization of simulated scenarios",
    colour = "True\npseudotime"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90")
  )

ggsave(
  filename = file.path(project_dir, "results", "supplementary", "Supplementary_Figure3_UMAP_scenarios.pdf"),
  plot = sup3,
  width = 12,
  height = 4
)

ggsave(
  filename = file.path(project_dir, "results", "supplementary", "Supplementary_Figure3_UMAP_scenarios.png"),
  plot = sup3,
  width = 12,
  height = 4,
  dpi = 300
)


############################################################
## 15. Save session information
############################################################

sink(file.path(project_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Analysis complete.")
message("Main figures saved in: ", file.path(project_dir, "results", "figures"))
message("Supplementary figures saved in: ", file.path(project_dir, "results", "supplementary"))
message("Metrics saved in: ", file.path(project_dir, "results", "metrics"))

