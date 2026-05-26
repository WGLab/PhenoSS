## ============================================================
##  PhenoSS - visualization pipeline
##  Embeddings : MDS | t-SNE | UMAP
##  Distances  : Resnik-only (alpha = 1) | PhenoSS-only (alpha = 0) | hybrid (0 < alpha < 1)
##  Overlay    : KMeans
##
##  PhenoSS DISTANCE:
##    Cross-patient contrastive distance using score x rank_frac
##    over all causal targets in the cohort. For patients i, j
##    with causal targets ci, cj:
##       w_p(t)   = -score[p, t] * rank_frac[p, t]
##       D[i, j] = -( w_i(cj) + w_j(ci) )
##    Rescaled to [0, 1].
##
##  OUTPUT:
##    Static ggplot2 figures only.
##    One PDF + PNG + TIFF per coord method.
##
##  N_CLUSTERS modes:
##    "true"  - k = unique true labels; clustering for stats only;
##              points colored and shaped by real labels exclusively.
##    "auto"  - algorithm decides
##    integer - force exactly k clusters
## ============================================================

library(vegan)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(Rtsne)
library(uwot)

# =============================================================================
# USER INPUTS
# =============================================================================
CASE_NAME   <- "casename here"
BASE_DIR    <- "modified here"
sim_file    <- file.path(BASE_DIR, paste0(CASE_NAME, "_phenoss_sim"))
scores_file <- file.path(BASE_DIR, "original_outputs.tsv")

TARGET_TYPE <- "gene"          # "gene" or "disease"

TSNE_PERPLEXITY  <- "auto"
TSNE_SEED        <- 42
UMAP_N_NEIGHBORS <- "auto"
UMAP_MIN_DIST    <- 0.1
UMAP_SEED        <- 42

COORD_METHODS <- c("MDS", "TSNE", "UMAP")
N_CLUSTERS    <- "true"

OUTLIER_IQR_THRESHOLD <- 0#3.0

PUB_WIDTH  <- 12.5
PUB_HEIGHT <- 6.5
PUB_DPI    <- 300

LABEL_TO_TARGET <- NULL
GROUP_ORDER     <- NULL

FLIP_DIM1_FOR <- c("TSNE", "UMAP")

VIZ_DIR <- file.path(BASE_DIR, "visualizations3")

# =============================================================================
# HELPERS
# =============================================================================
clean_label <- function(s) {
  s %>% as.character() %>% tolower() %>%
    str_replace_all("[_ ]", "-") %>% str_trim()
}

safe_rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

make_dot_sizes <- function(ranks, min_size = 1.5, max_size = 5) {
  ranks <- as.numeric(ranks)
  good  <- is.finite(ranks) & ranks > 0
  if (!any(good)) return(rep((min_size + max_size) / 2, length(ranks)))
  ranks[!good] <- max(ranks[good], na.rm = TRUE)
  log_inv <- log(1 / ranks)
  rng <- range(log_inv, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep((min_size + max_size) / 2, length(ranks)))
  }
  min_size + (max_size - min_size) * (log_inv - rng[1]) / (rng[2] - rng[1])
}

make_color_map <- function(groups) {
  base_cols <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
    "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"
  )
  cols <- if (length(groups) <= length(base_cols)) {
    base_cols[seq_along(groups)]
  } else {
    grDevices::rainbow(length(groups))
  }
  setNames(cols, groups)
}

LABEL_PCH <- c(16, 15, 17, 18, 8, 10, 7, 9, 11, 12, 13, 14, 3, 4, 6)
make_shape_map <- function(groups) {
  pchs <- LABEL_PCH[((seq_along(groups) - 1) %% length(LABEL_PCH)) + 1]
  setNames(pchs, groups)
}

trim_outliers_2d <- function(pts_2d, iqr_thresh = 3.0) {
  if (!is.finite(iqr_thresh) || iqr_thresh <= 0) return(rep(TRUE, nrow(pts_2d)))
  keep <- rep(TRUE, nrow(pts_2d))
  for (j in 1:2) {
    v   <- pts_2d[, j]
    q   <- quantile(v, c(0.25, 0.75), na.rm = TRUE)
    iqr <- q[2] - q[1]
    if (iqr == 0) next
    keep <- keep & (v >= q[1] - iqr_thresh * iqr) & (v <= q[2] + iqr_thresh * iqr)
  }
  keep
}

format_label_list <- function(labels) {
  labels <- as.character(labels)
  labels <- labels[nzchar(labels)]
  n <- length(labels)
  if (n == 0) return("")
  if (n == 1) return(labels)
  if (n == 2) return(paste(labels, collapse = " and "))
  paste0(paste(labels[-n], collapse = ", "), ", and ", labels[n])
}

prettify_label <- function(lbl) {
  if (exists("PRETTY_LABEL_MAP") && lbl %in% names(PRETTY_LABEL_MAP)) {
    return(unname(PRETTY_LABEL_MAP[[lbl]]))
  }
  tools::toTitleCase(gsub("-", " ", lbl))
}

make_cluster_color_map <- function(k) {
  pal <- c(
    "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
    "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
    "#CAB2D6", "#FFFF99", "#8DD3C7", "#BEBADA", "#80B1D3"
  )
  cols <- if (k <= length(pal)) pal[seq_len(k)] else grDevices::rainbow(k)
  setNames(cols, as.character(seq_len(k)))
}

compute_cluster_accuracy <- function(true_labels, pred_labels) {
  true_labels <- as.integer(as.factor(true_labels))
  pred_labels <- as.integer(as.factor(pred_labels))
  ct <- table(pred_labels, true_labels)
  sum(apply(ct, 1, max)) / length(true_labels)
}

compute_silhouette_score <- function(pts_2d, pred_labels) {
  if (!requireNamespace("cluster", quietly = TRUE)) return(NA_real_)
  ul <- sort(unique(pred_labels))
  if (length(ul) < 2) return(NA_real_)
  tryCatch({
    ss <- cluster::silhouette(pred_labels, dist(pts_2d))
    round(mean(ss[, 3]), 4)
  }, error = function(e) NA_real_)
}

resolve_k <- function(n_clusters_setting, n_true_labels) {
  if (identical(n_clusters_setting, "true")) return(n_true_labels)
  if (identical(n_clusters_setting, "auto")) return(NA_integer_)
  k <- suppressWarnings(as.integer(n_clusters_setting))
  if (is.na(k) || k < 2) stop("N_CLUSTERS must be 'true', 'auto', or integer >= 2")
  k
}

best_k_silhouette <- function(pts_2d, k_max) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Install 'cluster' for silhouette-based k selection.")
  }
  k_max  <- min(k_max, nrow(pts_2d) - 1)
  k_grid <- 2:k_max
  if (length(k_grid) == 0) return(2L)
  sils <- sapply(k_grid, function(k) {
    set.seed(42)
    km <- kmeans(pts_2d, centers = k, nstart = 10, iter.max = 200)
    ss <- cluster::silhouette(km$cluster, dist(pts_2d))
    mean(ss[, 3])
  })
  k_grid[which.max(sils)]
}

run_kmeans <- function(pts_2d, k_resolved, seed = 42) {
  if (is.na(k_resolved)) {
    k_resolved <- best_k_silhouette(pts_2d, k_max = floor(sqrt(nrow(pts_2d))))
    message(sprintf("    k-means auto k=%d (silhouette)", k_resolved))
  }
  set.seed(seed)
  kmeans(pts_2d, centers = k_resolved, nstart = 25, iter.max = 300)$cluster
}

build_meta <- function(x, ids, group_order = NULL) {
  m1 <- x[, c("id1", "label1")]
  colnames(m1) <- c("id", "group")
  m2 <- x[, c("id2", "label2")]
  colnames(m2) <- c("id", "group")
  
  meta <- rbind(m1, m2) %>%
    group_by(id) %>%
    summarize(
      group = {
        u <- unique(group)
        if (length(u) == 1) u else NA_character_
      },
      .groups = "drop"
    )
  
  meta <- data.frame(id = ids, stringsAsFactors = FALSE) %>%
    left_join(meta, by = "id")
  meta$group[is.na(meta$group)] <- "unknown"
  
  if (is.null(group_order)) {
    group_order <- sort(unique(meta$group))
  } else {
    extra <- setdiff(unique(meta$group), group_order)
    if (length(extra) > 0) group_order <- c(group_order, sort(extra))
  }
  
  meta$group <- factor(meta$group, levels = group_order)
  meta <- meta[match(ids, meta$id), , drop = FALSE]
  stopifnot(all(meta$id == ids))
  meta
}

resolve_target_column <- function(tt) {
  tt <- tolower(tt)
  if (tt == "gene") return("gene")
  if (tt == "disease") return("disease_mondo")
  stop("TARGET_TYPE must be 'gene' or 'disease'")
}

prepare_scores <- function(scores_raw, target_col, target_type) {
  needed <- c("patient_id", target_col, "score")
  miss   <- setdiff(needed, colnames(scores_raw))
  if (length(miss) > 0) stop("scores_file missing: ", paste(miss, collapse = ", "))
  
  scores_raw %>%
    {
      if (target_type == "gene") {
        distinct(., patient_id, .data[[target_col]], .keep_all = TRUE)
      } else {
        distinct(.)
      }
    } %>%
    group_by(patient_id) %>%
    mutate(
      rank_target = as.integer(min_rank(score)),
      n_targets   = n(),
      rank_frac   = ifelse(n_targets == 1, 1, 1 - (rank_target - 1) / (n_targets - 1))
    ) %>%
    ungroup() %>%
    arrange(patient_id, rank_target)
}

resolve_label_to_target <- function(groups, label_to_target, target_values) {
  tl <- setNames(target_values, clean_label(target_values))
  map <- if (is.null(label_to_target)) {
    tl
  } else {
    m <- label_to_target
    names(m) <- clean_label(names(m))
    m
  }
  
  mg <- setdiff(groups, names(map))
  if (length(mg) > 0) stop("Labels with no target mapping: ", paste(mg, collapse = ", "))
  
  mt <- setdiff(unique(unname(map[groups])), target_values)
  if (length(mt) > 0) stop("Mapped targets not in scores_file: ", paste(mt, collapse = ", "))
  map
}

build_patient_target_matrices <- function(phenoss, ids, target_col) {
  all_targets <- sort(unique(as.character(phenoss[[target_col]])))
  
  score_mat <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = length(all_targets),
    dimnames = list(ids, all_targets)
  )
  rank_mat <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = length(all_targets),
    dimnames = list(ids, all_targets)
  )
  
  for (i in seq_len(nrow(phenoss))) {
    pid <- phenoss$patient_id[i]
    tgt <- as.character(phenoss[[target_col]][i])
    if (pid %in% ids && tgt %in% all_targets) {
      score_mat[pid, tgt] <- phenoss$score[i]
      rank_mat[pid, tgt]  <- phenoss$rank_target[i]
    }
  }
  
  for (p in ids) {
    rr <- rank_mat[p, ]
    if (any(is.na(rr))) {
      worst <- if (all(is.na(rr))) length(all_targets) else max(rr, na.rm = TRUE)
      rank_mat[p, is.na(rr)] <- worst
    }
  }
  
  for (t in all_targets) {
    sc <- score_mat[, t]
    if (any(is.na(sc))) {
      mn <- if (all(is.na(sc))) 0 else mean(sc, na.rm = TRUE)
      score_mat[is.na(sc), t] <- mn
    }
  }
  
  rank_frac_mat <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = length(all_targets),
    dimnames = list(ids, all_targets)
  )
  for (p in ids) {
    rr <- rank_mat[p, ]
    n_cand <- length(rr)
    rank_frac_mat[p, ] <- if (n_cand == 1) 1 else 1 - (rr - 1) / (n_cand - 1)
  }
  
  list(score = score_mat, rank = rank_mat, rank_frac = rank_frac_mat, targets = all_targets)
}

build_contrastive_dist <- function(ptm, ids, patient_causal) {
  n <- length(ids)
  w_mat <- (-ptm$score) * ptm$rank_frac
  
  D <- matrix(0, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      ci <- patient_causal[ids[i]]
      cj <- patient_causal[ids[j]]
      if (is.na(ci) || is.na(cj) || !(ci %in% ptm$targets) || !(cj %in% ptm$targets)) {
        D[i, j] <- NA_real_
        next
      }
      w_ij <- w_mat[ids[i], cj]
      w_ji <- w_mat[ids[j], ci]
      D[i, j] <- -(w_ij + w_ji)
    }
  }
  
  if (any(is.na(D))) D[is.na(D)] <- max(D, na.rm = TRUE)
  diag(D) <- 0
  D <- (D + t(D)) / 2
  
  out <- safe_rescale01(D)
  rownames(out) <- ids
  colnames(out) <- ids
  out
}

embed_2d <- function(dist_mat, method, ids,
                     tsne_perplexity = "auto", tsne_seed = 42,
                     umap_n_neighbors = "auto", umap_min_dist = 0.1,
                     umap_seed = 42) {
  n <- nrow(dist_mat)
  method <- toupper(method)
  
  if (method == "MDS") {
    pts <- cmdscale(as.dist(dist_mat), k = 2, eig = TRUE)$points
  } else if (method == "TSNE") {
    perp <- if (identical(tsne_perplexity, "auto")) {
      max(1, min(floor(sqrt(n)), floor((n - 1) / 3)))
    } else {
      max(1, min(as.numeric(tsne_perplexity), floor((n - 1) / 3)))
    }
    message(sprintf("    t-SNE perplexity = %d  (n=%d)", perp, n))
    set.seed(tsne_seed)
    pts <- Rtsne(
      as.dist(dist_mat),
      dims = 2,
      perplexity = perp,
      is_distance = TRUE,
      check_duplicates = FALSE,
      verbose = FALSE
    )$Y
  } else if (method == "UMAP") {
    n_nb <- if (identical(umap_n_neighbors, "auto")) {
      max(5, min(15, floor(n / 10)))
    } else {
      min(as.integer(umap_n_neighbors), n - 1)
    }
    n_nb <- min(n_nb, n - 1)
    message(sprintf("    UMAP n_neighbors = %d  (n=%d)", n_nb, n))
    dist_m  <- as.matrix(dist_mat)
    nn_idx  <- t(apply(dist_m, 1, function(row) order(row)[2:(n_nb + 1)]))
    nn_dist <- t(apply(dist_m, 1, function(row) sort(row)[2:(n_nb + 1)]))
    set.seed(umap_seed)
    pts <- uwot::umap(
      X = NULL,
      nn_method = list(idx = nn_idx, dist = nn_dist),
      n_components = 2,
      min_dist = umap_min_dist,
      verbose = FALSE
    )
  } else {
    stop("Unknown method: ", method)
  }
  
  rownames(pts) <- ids
  pts
}

distance_label <- function(key) {
  if (identical(key, "resnik_only")) return("Resnik only (alpha = 1.0)")
  if (identical(key, "phenoss_only")) return("PhenoSS only (alpha = 0.0)")
  key
}

compute_kmeans_stats <- function(pts_2d, k_resolved, true_int, has_true) {
  km_labels <- tryCatch(run_kmeans(pts_2d, k_resolved), error = function(e) NULL)
  list(
    km_labels = km_labels,
    km_acc = if (has_true && !is.null(km_labels)) compute_cluster_accuracy(true_int, km_labels) else NULL,
    km_sil = if (!is.null(km_labels)) compute_silhouette_score(pts_2d, km_labels) else NA_real_
  )
}

make_static_subplot <- function(pts_2d, groups, dot_sizes, present_groups,
                                color_map, shape_map,
                                R2_val, p_val, dist_label,
                                km_labels = NULL, km_acc = NULL, km_sil = NULL,
                                n_clusters_setting = "true",
                                outlier_iqr = 3.0,
                                show_legend = FALSE) {
  df <- data.frame(
    x = pts_2d[, 1],
    y = pts_2d[, 2],
    group = factor(groups, levels = present_groups),
    size = dot_sizes,
    stringsAsFactors = FALSE
  )
  
  keep <- trim_outliers_2d(pts_2d, iqr_thresh = outlier_iqr)
  df <- df[keep, ]
  
  group_counts <- table(df$group)
  df$group <- factor(df$group, levels = names(sort(group_counts, decreasing = TRUE)))
  df <- df[order(df$group), ]
  
  xpad <- diff(range(df$x)) * 0.04
  ypad <- diff(range(df$y)) * 0.04
  xlim_plot <- range(df$x) + c(-xpad, xpad)
  ylim_plot <- range(df$y) + c(-ypad, ypad)
  
  p <- ggplot(df, aes(x = x, y = y, color = group, shape = group, size = size))
  
  if (!is.null(km_labels) && !identical(n_clusters_setting, "true")) {
    cl_df <- data.frame(
      x = pts_2d[keep, 1],
      y = pts_2d[keep, 2],
      cluster = factor(km_labels[keep])
    )
    cluster_colors <- make_cluster_color_map(max(as.integer(levels(cl_df$cluster))))
    for (cl in levels(cl_df$cluster)) {
      sub <- cl_df[cl_df$cluster == cl, ]
      if (nrow(sub) < 3) next
      hull_idx <- chull(sub$x, sub$y)
      hull_df <- sub[c(hull_idx, hull_idx[1]), ]
      p <- p +
        geom_polygon(
          data = hull_df,
          aes(x = x, y = y),
          inherit.aes = FALSE,
          fill = cluster_colors[cl],
          alpha = 0.12,
          color = NA
        )
    }
  }
  
  p <- p +
    geom_point(alpha = 0.88, stroke = 0.3) +
    scale_color_manual(values = color_map, name = "Group", drop = FALSE) +
    scale_shape_manual(values = shape_map, name = "Group", drop = FALSE) +
    scale_size_identity() +
    coord_cartesian(xlim = xlim_plot, ylim = ylim_plot, expand = FALSE) +
    labs(x = "Dim 1", y = "Dim 2", title = dist_label) +
    theme_classic(base_size = 11) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.background = element_rect(fill = "white", color = "grey60", linewidth = 0.3),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = grid::unit(0.55, "cm"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  
  stat_lines <- c(
    sprintf("PERMANOVA R2=%.2f  p=%s", R2_val, signif(p_val, 3)),
    if (!is.null(km_acc) && !is.na(km_acc)) sprintf("Cluster acc = %.4f", km_acc),
    if (!is.null(km_sil) && !is.na(km_sil)) sprintf("Silhouette = %.4f", km_sil)
  )
  stat_label <- paste(stat_lines, collapse = "\n")
  
  p +
    annotate(
      "label",
      x = xlim_plot[1] + 0.02 * diff(xlim_plot),
      y = ylim_plot[1] + 0.02 * diff(ylim_plot),
      label = stat_label,
      hjust = 0,
      vjust = 0,
      size = 3.5,
      lineheight = 1.5,
      label.padding = grid::unit(0.45, "lines"),
      label.r = grid::unit(0.12, "lines"),
      fill = "white",
      color = "black",
      fontface = "plain",
      alpha = 0.93
    )
}

# =============================================================================
# STEP 1 - SETUP
# =============================================================================
TARGET_COL <- resolve_target_column(TARGET_TYPE)
message("TARGET_TYPE=", TARGET_TYPE, "  column=", TARGET_COL)
if (!dir.exists(VIZ_DIR)) {
  dir.create(VIZ_DIR, recursive = TRUE)
  message("Created: ", VIZ_DIR)
}

# =============================================================================
# STEP 2 - RESNIK SIMILARITY FILE
# =============================================================================
x <- read.table(sim_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(x) <- c("id1", "id2", "sim", "label1", "label2")
x <- x[!is.na(x$sim), ]
x$id1 <- trimws(x$id1)
x$id2 <- trimws(x$id2)

raw_pretty <- unique(c(trimws(as.character(x$label1)), trimws(as.character(x$label2))))
PRETTY_LABEL_MAP <- setNames(raw_pretty, clean_label(raw_pretty))

x$label1 <- clean_label(x$label1)
x$label2 <- clean_label(x$label2)
x$dist   <- pmin(1 / x$sim, 1)
stopifnot(all(x$dist > 0), all(x$dist <= 1), all(is.finite(x$dist)))

# =============================================================================
# STEP 3 - PHENOSS SCORES
# =============================================================================
scores_raw <- read.table(
  scores_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)
phenoss <- prepare_scores(scores_raw, TARGET_COL, TARGET_TYPE)
target_values <- sort(unique(as.character(phenoss[[TARGET_COL]])))
message("Targets: ", paste(target_values, collapse = ", "))

# =============================================================================
# STEP 4 - IDS + META
# =============================================================================
ids <- sort(unique(c(x$id1, x$id2)))
n <- length(ids)
raw_groups <- sort(unique(c(x$label1, x$label2)))
has_true_labels_global <- length(raw_groups) > 0 && !all(raw_groups %in% c("", "na", "unknown"))
n_true_labels <- if (has_true_labels_global) length(raw_groups) else 0L

label_to_target_map <- resolve_label_to_target(raw_groups, LABEL_TO_TARGET, target_values)
if (is.null(GROUP_ORDER)) GROUP_ORDER <- sort(raw_groups)

meta <- build_meta(x, ids, group_order = GROUP_ORDER)
present_groups <- levels(droplevels(meta$group))
true_int <- as.integer(meta$group)

message("Groups (", n_true_labels, "): ", paste(present_groups, collapse = ", "))

k_resolved <- resolve_k(N_CLUSTERS, n_true_labels)
message(sprintf(
  "N_CLUSTERS='%s' -> k_resolved=%s",
  N_CLUSTERS,
  if (is.na(k_resolved)) "auto" else k_resolved
))

# =============================================================================
# STEP 5 - COLORS + SHAPES
# =============================================================================
color_map <- make_color_map(present_groups)
shape_map <- make_shape_map(present_groups)

# =============================================================================
# STEP 6 - RESNIK DISTANCE MATRIX
# =============================================================================
dist_mat_pheno <- matrix(1, n, n, dimnames = list(ids, ids))
diag(dist_mat_pheno) <- 0
for (i in seq_len(nrow(x))) {
  dist_mat_pheno[x$id1[i], x$id2[i]] <- x$dist[i]
  dist_mat_pheno[x$id2[i], x$id1[i]] <- x$dist[i]
}

# =============================================================================
# STEP 7 - CONTRASTIVE PHENOSS DISTANCE
# =============================================================================
ptm <- build_patient_target_matrices(phenoss, ids, TARGET_COL)
message(sprintf(
  "Patient x target matrix: %d patients x %d targets",
  nrow(ptm$score),
  ncol(ptm$score)
))

patient_causal <- unname(label_to_target_map[as.character(meta$group)])
names(patient_causal) <- ids

dist_mat_target <- build_contrastive_dist(ptm, ids, patient_causal)
message("Contrastive PhenoSS distance built.")

chosen_stats <- data.frame(
  patient_id = ids,
  target_label = as.character(meta$group),
  target_value = unname(patient_causal[ids]),
  score = ptm$score[cbind(ids, unname(patient_causal[ids]))],
  rank_target = ptm$rank[cbind(ids, unname(patient_causal[ids]))],
  stringsAsFactors = FALSE
)
chosen_stats$rank_frac <- ptm$rank_frac[cbind(
  chosen_stats$patient_id,
  chosen_stats$target_value
)]

# =============================================================================
# STEP 8 - DOT SIZES
# =============================================================================
proband_rank <- setNames(chosen_stats$rank_target, chosen_stats$patient_id)[ids]
dot_sizes_gg <- make_dot_sizes(proband_rank, min_size = 1.5, max_size = 5)

# =============================================================================
# STEP 9 - PRE-COMPUTE PERMANOVA + DISTANCE SETS
# =============================================================================
message("\n=== Pre-computing PERMANOVA stats ===")

distance_sets <- list(
  resnik_only = dist_mat_pheno,
  phenoss_only = dist_mat_target
)

perm_stats <- lapply(names(distance_sets), function(key) {
  message("  ", key, " ...")
  mat <- distance_sets[[key]]
  diag(mat) <- 0
  mat <- (mat + t(mat)) / 2
  pa <- adonis2(as.dist(mat) ~ group, data = meta, permutations = 999)
  list(
    mat = mat,
    R2 = round(pa$R2[1], 2),
    p = pa$`Pr(>F)`[1]
  )
})
names(perm_stats) <- names(distance_sets)

# =============================================================================
# STEP 10 - STATIC PUBLICATION FIGURES
# =============================================================================
message("\n=== Generating static publication figures ===")

pretty_groups <- vapply(present_groups, prettify_label, character(1))
cohort_str <- format_label_list(pretty_groups)

for (cm in COORD_METHODS) {
  message(sprintf("  Coord method: %s", cm))
  
  embed_cache <- lapply(names(perm_stats), function(key) {
    message(sprintf("    Embedding: %s | %s", key, cm))
    pts_2d <- embed_2d(
      perm_stats[[key]]$mat,
      cm,
      ids,
      tsne_perplexity = TSNE_PERPLEXITY,
      tsne_seed = TSNE_SEED,
      umap_n_neighbors = UMAP_N_NEIGHBORS,
      umap_min_dist = UMAP_MIN_DIST,
      umap_seed = UMAP_SEED
    )
    
    if (cm %in% FLIP_DIM1_FOR) {
      pts_2d[, 1] <- -pts_2d[, 1]
    }
    
    has_true <- n_true_labels > 0
    km_stats <- compute_kmeans_stats(pts_2d, k_resolved, true_int, has_true)
    
    list(
      pts_2d = pts_2d,
      R2_val = perm_stats[[key]]$R2,
      p_val = perm_stats[[key]]$p,
      dist_label = distance_label(key),
      km_labels = km_stats$km_labels,
      km_acc = km_stats$km_acc,
      km_sil = km_stats$km_sil
    )
  })
  names(embed_cache) <- names(perm_stats)
  
  subplots <- list(
    make_static_subplot(
      pts_2d = embed_cache[["resnik_only"]]$pts_2d,
      groups = as.character(meta$group),
      dot_sizes = dot_sizes_gg,
      present_groups = present_groups,
      color_map = color_map,
      shape_map = shape_map,
      R2_val = embed_cache[["resnik_only"]]$R2_val,
      p_val = embed_cache[["resnik_only"]]$p_val,
      dist_label = embed_cache[["resnik_only"]]$dist_label,
      km_labels = embed_cache[["resnik_only"]]$km_labels,
      km_acc = embed_cache[["resnik_only"]]$km_acc,
      km_sil = embed_cache[["resnik_only"]]$km_sil,
      n_clusters_setting = N_CLUSTERS,
      outlier_iqr = OUTLIER_IQR_THRESHOLD,
      show_legend = FALSE
    ),
    make_static_subplot(
      pts_2d = embed_cache[["phenoss_only"]]$pts_2d,
      groups = as.character(meta$group),
      dot_sizes = dot_sizes_gg,
      present_groups = present_groups,
      color_map = color_map,
      shape_map = shape_map,
      R2_val = embed_cache[["phenoss_only"]]$R2_val,
      p_val = embed_cache[["phenoss_only"]]$p_val,
      dist_label = embed_cache[["phenoss_only"]]$dist_label,
      km_labels = embed_cache[["phenoss_only"]]$km_labels,
      km_acc = embed_cache[["phenoss_only"]]$km_acc,
      km_sil = embed_cache[["phenoss_only"]]$km_sil,
      n_clusters_setting = N_CLUSTERS,
      outlier_iqr = OUTLIER_IQR_THRESHOLD,
      show_legend = TRUE
    )
  )
  
  combined <- wrap_plots(subplots, nrow = 1, guides = "collect") +
    plot_annotation(
      title = sprintf(
        "Patient Clustering of %s Syndromes  |  %s  |  KMeans",
        cohort_str, cm
      ),
      theme = theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
  
  base <- sprintf("%s_%s_resnik_vs_phenoss_kmeans", CASE_NAME, cm)
  pdf_path  <- file.path(VIZ_DIR, paste0(base, ".pdf"))
  png_path  <- file.path(VIZ_DIR, paste0(base, ".png"))
  tiff_path <- file.path(VIZ_DIR, paste0(base, ".tiff"))
  
  ggsave(
    pdf_path, combined,
    width = PUB_WIDTH, height = PUB_HEIGHT,
    units = "in", device = cairo_pdf
  )
  message("    Saved PDF  -> ", pdf_path)
  
  ggsave(
    png_path, combined,
    width = PUB_WIDTH, height = PUB_HEIGHT,
    units = "in", dpi = PUB_DPI
  )
  message("    Saved PNG  -> ", png_path)
  
  ggsave(
    tiff_path, combined,
    width = PUB_WIDTH, height = PUB_HEIGHT,
    units = "in", dpi = PUB_DPI, device = "tiff", compression = "lzw"
  )
  message("    Saved TIFF -> ", tiff_path)
}

message("\n=== Final figures done ===")