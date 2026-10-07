## ============================================================
##  PhenoSS — Patient Clustering & Visualization Pipeline
##  Embeddings : MDS | t-SNE | UMAP
##  Distances  : Resnik-only | PhenoSS-only | Combined (alpha blend)
##  Overlays   : Voronoi | Density | HDBSCAN | KMeans
##
##  QUICK-START — three things you must set before running:
##    1. sim_file    — see "INPUT FILE FORMATS" below
##    2. scores_file — the run_phenoss.py output TSV
##    3. TARGET_TYPE — "gene" or "disease"
##
##  Everything else in USER INPUTS has a working default.
##
## ============================================================
##
##  INPUT FILE FORMATS
##  ------------------
##
##  sim_file  (tab-separated, NO header row):
##    Column 1: patient_id_1       (string)
##    Column 2: patient_id_2       (string)
##    Column 3: similarity_score   (numeric, > 0)
##    Column 4: label_1            (group/disease label for patient 1)
##    Column 5: label_2            (group/disease label for patient 2)
##
##    IMPORTANT: this is NOT the raw phenoss_similarity.tsv produced by
##    similarity_score.py (which has HPO terms in columns 3-4, not labels).
##    You need to join your patient->label mapping onto the similarity
##    output first. Example in Python:
##
##      import pandas as pd
##      CASE_NAME = 'PUT IT HERE'
##      BASE_DIR = 'PUT IT HERE'
##      sim   = pd.read_csv("phenoss_similarity.tsv", sep="\t")
##      labels = pd.read_csv("patient_labels.tsv", sep="\t")
##      # labels has columns: patient_id, label
##      sim = sim.merge(labels.rename(columns={"patient_id":"pat1","label":"label1"}), on="pat1")
##      sim = sim.merge(labels.rename(columns={"patient_id":"pat2","label":"label2"}), on="pat2")
##      sim[["pat1","pat2","similarity","label1","label2"]].to_csv(
##          f"{BASE_DIR}/{CASE_NAME}_phenoss_sim"), sep="\t", index=False, header=False)
##
##    If your labels ARE already the disease/gene IDs in the scores file,
##    set LABEL_TO_TARGET <- NULL below (labels = targets directly).
##
##  scores_file  (tab-separated, WITH header row):
##    Required columns (from run_phenoss.py output):
##      patient_id, score, and either:
##        disease_mondo  (when TARGET_TYPE = "disease")
##        gene           (when TARGET_TYPE = "gene")
##    Additional columns (rank, disease_id, etc.) are ignored.
##
## ============================================================
##
##  ALPHA — what it means
##  ---------------------
##  The combined distance is:
##    D = alpha * D_resnik + (1 - alpha) * D_phenoss
##
##    alpha = 1.0  ->  Resnik semantic similarity only  (HPO structure)
##    alpha = 0.0  ->  PhenoSS only  (recommended; uses gene/disease scores)
##    alpha = 0.5  ->  equal blend of both
##
##  The script runs PERMANOVA for every alpha in ALPHA_VALUES and saves
##  separate panels for each alpha listed in PUB_ALPHAS. You do NOT need
##  to pick a single alpha — all are computed and shown side-by-side.
##  To visualize only one alpha, set PUB_ALPHAS <- c(your_alpha).
##
## ============================================================

## ---- dependencies ---------------------------------------------------
library(vegan)        # adonis2 (PERMANOVA)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)    # combine subplots
library(Rtsne)
library(uwot)
library(dbscan)

## ---- plotly (needed only for interactive HTML loop) -----------------
# library(plotly)
# library(htmlwidgets)

# =============================================================================
# USER INPUTS — REQUIRED
# =============================================================================

## ---- File paths ---------------------------------------------------------
CASE_NAME   <- "casename here"
BASE_DIR    <- "modified here"

sim_file    <- file.path(BASE_DIR, paste0(CASE_NAME, "_phenoss_sim"))
# ^ Tab-separated, no header. Format: id1  id2  similarity  label1  label2
# See INPUT FILE FORMATS at the top of this script.

scores_file <- file.path(BASE_DIR, "original_outputs.tsv")
# ^ run_phenoss.py output. Must have columns: patient_id, score, and
#   disease_mondo (if TARGET_TYPE="disease") or gene (if TARGET_TYPE="gene").

## ---- Target type -------------------------------------------------------
TARGET_TYPE <- "disease"   # "gene" or "disease"
# "disease" -> matches the disease_mondo column in scores_file
# "gene"    -> matches the gene column in scores_file

## ---- Label -> Target mapping -------------------------------------------
# This tells the script how to connect sim_file labels to scores_file targets.
#
# NULL means the labels in sim_file ARE the gene/MONDO IDs in scores_file.
# Use NULL when, for example, label1/label2 are already "MONDO:0007947" etc.
#
# Otherwise, provide a named character vector mapping each label to its target:
#   LABEL_TO_TARGET <- c(
#     "taf1-associated"   = "TAF1",           # for TARGET_TYPE="gene"
#     "cornelia-de-lange" = "MONDO:0007947"   # for TARGET_TYPE="disease"
#   )
# Labels are matched case-insensitively after replacing spaces/underscores
# with hyphens (the clean_label() function). Target values must match exactly
# what appears in the scores_file target column.
LABEL_TO_TARGET <- NULL

# =============================================================================
# USER INPUTS — ALPHA (distance blend)
# =============================================================================

## ---- Alpha grid (all values that get computed and plotted) --------------
# "resnik_only" is a special token meaning alpha=1.0 (no PhenoSS component).
# Add or remove numeric alpha values as needed.
ALPHA_VALUES <- c("resnik_only", 0.0, 0.5)
# Common choices:
#   c("resnik_only", 0.0)         -> just compare Resnik vs PhenoSS-only
#   c("resnik_only", 0.0, 0.5)    -> add the 50/50 blend
#   c("resnik_only", 0.0, 0.4, 0.5, 0.6)  -> sweep around 0.5

## ---- Publication panels (subset of ALPHA_VALUES shown in saved figures) -
# Left -> right order. Each entry must appear in ALPHA_VALUES above,
# or be "resnik_only". Typically 3 panels fit well at PUB_WIDTH=18in.
PUB_ALPHAS <- c("resnik_only", 0.5, 0.0)

# =============================================================================
# USER INPUTS — VISUALIZATION
# =============================================================================

## ---- Embedding methods --------------------------------------------------
COORD_METHODS <- c("MDS", "TSNE", "UMAP")   # remove any you don't want

## ---- Overlay style ------------------------------------------------------
# One figure is saved per (coord method, overlay). Choose any subset:
#   "voronoi"  — Voronoi region fill by group centroid
#   "density"  — 2D density contours per group
#   "hdbscan"  — HDBSCAN cluster convex hulls
#   "kmeans"   — k-means cluster convex hulls
PUB_OVERLAYS <- c("kmeans", "voronoi")

## ---- Clustering K -------------------------------------------------------
# Controls k for HDBSCAN and k-means overlays (not for point coloring,
# which always uses the true labels from sim_file).
#   "true"    -> k = number of distinct labels in sim_file (recommended)
#   "auto"    -> algorithm picks k via silhouette score
#   integer   -> force exactly this many clusters (e.g. N_CLUSTERS <- 3)
N_CLUSTERS <- "true"

## ---- t-SNE settings -----------------------------------------------------
TSNE_PERPLEXITY <- "auto"   # "auto" = sqrt(n), or set a numeric value
TSNE_SEED       <- 42

## ---- UMAP settings ------------------------------------------------------
UMAP_N_NEIGHBORS <- "auto"   # "auto" = clamp(n/10, 5, 15), or set numeric
UMAP_MIN_DIST    <- 0.1
UMAP_SEED        <- 42

## ---- Figure dimensions --------------------------------------------------
PUB_WIDTH  <- 18    # inches
PUB_HEIGHT <- 6.5
PUB_DPI    <- 300

## ---- Outlier trimming ---------------------------------------------------
# Points > N IQR-widths from Q1/Q3 in either embedding dimension are hidden.
# Set Inf to disable trimming entirely.
OUTLIER_IQR_THRESHOLD <- 3.0

# =============================================================================
# USER INPUTS — ADVANCED (safe to leave as-is)
# =============================================================================

## ---- Group display order ------------------------------------------------
# NULL = alphabetical. Set to a character vector to control legend/color order.
# Example: GROUP_ORDER <- c("taf1-associated", "brd4", "nipbl", "smc3")
GROUP_ORDER <- NULL

## ---- Voronoi region bias ------------------------------------------------
# NULL = equal bias for all groups. Named numeric vector to expand/shrink
# regions for specific groups, e.g.: GROUP_BIAS <- c("smc3" = 1.5)
GROUP_BIAS <- NULL

## ---- Voronoi grid resolution --------------------------------------------
GRID_RES <- 800   # higher = sharper boundaries, slower
PAD_FRAC <- 0.15  # fraction of axis range added as padding around points

## ---- Output directory ---------------------------------------------------
VIZ_DIR <- file.path(BASE_DIR, "visualizations")

# =============================================================================
# HELPERS
# =============================================================================
clean_label <- function(s) {
  s %>% as.character() %>% tolower() %>%
    str_replace_all("[_ ]", "-") %>% str_trim()
}

safe_rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0)
    return(matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

make_dot_sizes <- function(ranks, min_size = 1.5, max_size = 5) {
  ranks <- as.numeric(ranks)
  good  <- is.finite(ranks) & ranks > 0
  if (!any(good)) return(rep((min_size + max_size) / 2, length(ranks)))
  ranks[!good] <- max(ranks[good], na.rm = TRUE)
  log_inv <- log(1 / ranks)
  rng <- range(log_inv, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0)
    return(rep((min_size + max_size) / 2, length(ranks)))
  min_size + (max_size - min_size) * (log_inv - rng[1]) / (rng[2] - rng[1])
}

make_color_map <- function(groups) {
  base_cols <- c(
    "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
    "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
    "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
    "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E"
  )
  cols <- if (length(groups) <= length(base_cols))
    base_cols[seq_along(groups)] else grDevices::rainbow(length(groups))
  setNames(cols, groups)
}

LABEL_PCH <- c(16, 15, 17, 18, 8, 10, 7, 9, 11, 12, 13, 14, 3, 4, 6)
make_shape_map <- function(groups) {
  pchs <- LABEL_PCH[((seq_along(groups) - 1) %% length(LABEL_PCH)) + 1]
  setNames(pchs, groups)
}

LABEL_SYMBOLS <- c(
  "circle","square","diamond","cross","x",
  "triangle-up","triangle-down","triangle-left",
  "triangle-right","pentagon","hexagon","star",
  "hourglass","bowtie","asterisk"
)
make_symbol_map <- function(groups) {
  syms <- LABEL_SYMBOLS[((seq_along(groups) - 1) %% length(LABEL_SYMBOLS)) + 1]
  setNames(syms, groups)
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

# ---------- Dynamic title helpers ----------------------------------------
format_label_list <- function(labels) {
  labels <- as.character(labels); labels <- labels[nzchar(labels)]
  n <- length(labels)
  if (n == 0) return("")
  if (n == 1) return(labels)
  if (n == 2) return(paste(labels, collapse = " and "))
  paste0(paste(labels[-n], collapse = ", "), ", and ", labels[n])
}

prettify_label <- function(lbl) {
  # Uses PRETTY_LABEL_MAP built automatically from raw sim_file labels (STEP 2).
  # Fallback: title-case the cleaned label.
  if (exists("PRETTY_LABEL_MAP") && lbl %in% names(PRETTY_LABEL_MAP))
    return(unname(PRETTY_LABEL_MAP[[lbl]]))
  tools::toTitleCase(gsub("-", " ", lbl))
}

# ---------- Voronoi helpers ----------------------------------------------
compute_voronoi_objects <- function(pts_2d, groups, present_groups,
                                    grid_res = 800, pad_frac = 0.15,
                                    group_bias = NULL) {
  df        <- data.frame(Dim1 = pts_2d[,1], Dim2 = pts_2d[,2], group = groups)
  centroids <- df %>% group_by(group) %>%
    summarize(cx = mean(Dim1), cy = mean(Dim2), .groups = "drop")
  xr <- range(pts_2d[,1]); yr <- range(pts_2d[,2])
  xp <- diff(xr) * pad_frac; yp <- diff(yr) * pad_frac
  xlim <- c(xr[1] - xp, xr[2] + xp); ylim <- c(yr[1] - yp, yr[2] + yp)
  x_seq <- seq(xlim[1], xlim[2], length.out = grid_res)
  y_seq <- seq(ylim[1], ylim[2], length.out = grid_res)
  gdf   <- expand.grid(Dim1 = x_seq, Dim2 = y_seq)
  if (is.null(group_bias)) {
    group_bias <- rep(1.0, length(present_groups)); names(group_bias) <- present_groups
  } else {
    tmp <- rep(1.0, length(present_groups)); names(tmp) <- present_groups
    tmp[names(group_bias)] <- group_bias; group_bias <- tmp[present_groups]
  }
  d2c <- sapply(seq_len(nrow(centroids)), function(k) {
    gi <- match(as.character(centroids$group[k]), present_groups)
    sqrt((gdf$Dim1 - centroids$cx[k])^2 + (gdf$Dim2 - centroids$cy[k])^2) / group_bias[gi]
  })
  nearest   <- apply(d2c, 1, which.min)
  label_int <- match(as.character(centroids$group[nearest]), present_groups)
  z_mat     <- matrix(label_int, nrow = grid_res, ncol = grid_res)
  pt_d2c    <- sapply(seq_len(nrow(centroids)), function(k)
    sqrt((df$Dim1 - centroids$cx[k])^2 + (df$Dim2 - centroids$cy[k])^2))
  pt_assigned <- as.character(centroids$group[apply(pt_d2c, 1, which.min)])
  region_acc  <- mean(pt_assigned == as.character(groups))
  list(centroids = centroids, x_seq = x_seq, y_seq = y_seq,
       z_mat = z_mat, xlim = xlim, ylim = ylim, region_acc = region_acc)
}

get_boundary_shapes <- function(x_seq, y_seq, z_mat) {
  levs <- seq_len(length(unique(as.vector(z_mat))) - 1) + 0.5
  cls  <- contourLines(x = x_seq, y = y_seq, z = z_mat, levels = levs)
  lapply(cls, function(cl) {
    path <- paste0("M ", cl$x[1], ",", cl$y[1], " ",
                   paste(sprintf("L %s,%s", cl$x[-1], cl$y[-1]), collapse = " "))
    list(type = "path", layer = "above", path = path,
         line = list(color = "rgba(0,0,0,0.7)", width = 2))
  })
}

make_voronoi_rects <- function(x_seq, y_seq, z_mat, present_groups,
                               group_colors, alpha = 0.20, ds = 6) {
  idx_x <- seq(1, length(x_seq), by = ds)
  idx_y <- seq(1, length(y_seq), by = ds)
  z_ds  <- z_mat[idx_x, idx_y, drop = FALSE]
  xs    <- x_seq[idx_x]; ys <- y_seq[idx_y]
  dx    <- if (length(xs) > 1) diff(xs)[1] else diff(x_seq)[1] * ds
  dy    <- if (length(ys) > 1) diff(ys)[1] else diff(y_seq)[1] * ds
  shapes <- list()
  for (i in seq_len(nrow(z_ds))) for (j in seq_len(ncol(z_ds))) {
    grp  <- present_groups[z_ds[i, j]]
    hex  <- group_colors[grp]; rgbv <- col2rgb(hex)
    fill <- sprintf("rgba(%d,%d,%d,%.2f)", rgbv[1], rgbv[2], rgbv[3], alpha)
    shapes[[length(shapes) + 1]] <- list(
      type = "rect", layer = "below",
      x0 = xs[i] - dx/2, x1 = xs[i] + dx/2,
      y0 = ys[j] - dy/2, y1 = ys[j] + dy/2,
      fillcolor = fill, line = list(width = 0))
  }
  shapes
}

# ---------- Clustering helpers -------------------------------------------
make_cluster_color_map <- function(k) {
  pal <- c("#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A",
           "#B15928","#A6CEE3","#B2DF8A","#FB9A99","#FDBF6F",
           "#CAB2D6","#FFFF99","#8DD3C7","#BEBADA","#80B1D3")
  cols <- if (k <= length(pal)) pal[seq_len(k)] else grDevices::rainbow(k)
  setNames(cols, as.character(seq_len(k)))
}

compute_cluster_accuracy <- function(true_labels, pred_labels) {
  true_labels <- as.integer(as.factor(true_labels))
  pred_labels <- as.integer(as.factor(pred_labels))
  ct  <- table(pred_labels, true_labels)
  sum(apply(ct, 1, max)) / length(true_labels)
}

compute_silhouette_score <- function(pts_2d, pred_labels) {
  if (!requireNamespace("cluster", quietly = TRUE)) return(NA_real_)
  ul <- sort(unique(pred_labels[pred_labels != 0]))
  if (length(ul) < 2) return(NA_real_)
  valid <- pred_labels != 0
  if (sum(valid) < 3) return(NA_real_)
  tryCatch({
    ss <- cluster::silhouette(pred_labels[valid], dist(pts_2d[valid, , drop = FALSE]))
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
  if (!requireNamespace("cluster", quietly = TRUE))
    stop("Install 'cluster' for silhouette-based k selection.")
  k_max  <- min(k_max, nrow(pts_2d) - 1)
  k_grid <- 2:k_max
  if (length(k_grid) == 0) return(2L)
  sils <- sapply(k_grid, function(k) {
    set.seed(42)
    km  <- kmeans(pts_2d, centers = k, nstart = 10, iter.max = 200)
    ss  <- cluster::silhouette(km$cluster, dist(pts_2d))
    mean(ss[, 3])
  })
  k_grid[which.max(sils)]
}

run_hdbscan <- function(pts_2d, k_resolved) {
  minPts <- max(3, floor(log(nrow(pts_2d))))
  labels <- dbscan::hdbscan(pts_2d, minPts = minPts)$cluster
  if (!is.na(k_resolved)) {
    non_noise <- labels[labels != 0]
    if (length(unique(non_noise)) > k_resolved) {
      cl_sizes <- sort(table(non_noise))
      while (length(unique(labels[labels != 0])) > k_resolved) {
        smallest <- as.integer(names(cl_sizes)[1]); cl_sizes <- cl_sizes[-1]
        labels[labels == smallest] <- 0
      }
    }
    if (any(labels == 0)) {
      active_cls <- sort(unique(labels[labels != 0]))
      if (length(active_cls) == 0) { labels[] <- 1L } else {
        centroids <- do.call(rbind, lapply(active_cls, function(cl)
          colMeans(pts_2d[labels == cl, , drop = FALSE])))
        for (i in which(labels == 0)) {
          dists <- rowSums((centroids - matrix(pts_2d[i,], nrow(centroids), 2, byrow = TRUE))^2)
          labels[i] <- active_cls[which.min(dists)]
        }
      }
    }
    labels <- as.integer(as.factor(labels))
  }
  labels
}

run_kmeans <- function(pts_2d, k_resolved, seed = 42) {
  if (is.na(k_resolved)) {
    k_resolved <- best_k_silhouette(pts_2d, k_max = floor(sqrt(nrow(pts_2d))))
    message(sprintf("    k-means auto k=%d (silhouette)", k_resolved))
  }
  set.seed(seed)
  kmeans(pts_2d, centers = k_resolved, nstart = 25, iter.max = 300)$cluster
}

# ---------- Meta / scores helpers ----------------------------------------
build_meta <- function(x, ids, group_order = NULL) {
  m1 <- x[, c("id1","label1")]; colnames(m1) <- c("id","group")
  m2 <- x[, c("id2","label2")]; colnames(m2) <- c("id","group")
  meta <- rbind(m1, m2) %>%
    group_by(id) %>%
    summarize(group = { u <- unique(group); if (length(u) == 1) u else NA_character_ },
              .groups = "drop")
  meta <- data.frame(id = ids, stringsAsFactors = FALSE) %>% left_join(meta, by = "id")
  meta$group[is.na(meta$group)] <- "unknown"
  if (is.null(group_order)) group_order <- sort(unique(meta$group))
  else {
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
  if (tt == "gene")    return("gene")
  if (tt == "disease") return("disease_mondo")
  stop("TARGET_TYPE must be 'gene' or 'disease'")
}

prepare_scores <- function(scores_raw, target_col, target_type) {
  needed <- c("patient_id", target_col, "score")
  miss   <- setdiff(needed, colnames(scores_raw))
  if (length(miss) > 0) stop(
    "scores_file is missing required columns: ", paste(miss, collapse = ", "), "\n",
    "  Found columns: ", paste(colnames(scores_raw), collapse = ", "), "\n",
    "  Check that TARGET_TYPE matches your scores_file format.\n",
    "  TARGET_TYPE='disease' requires column 'disease_mondo'.\n",
    "  TARGET_TYPE='gene'    requires column 'gene'.")
  scores_raw %>%
    { if (target_type == "gene") distinct(., patient_id, .data[[target_col]], .keep_all = TRUE)
      else distinct(.) } %>%
    group_by(patient_id) %>%
    mutate(rank_target = as.integer(min_rank(score)),
           n_targets   = n(),
           rank_frac   = ifelse(n_targets == 1, 1, 1 - (rank_target - 1) / (n_targets - 1))) %>%
    ungroup() %>%
    arrange(patient_id, rank_target)
}

resolve_label_to_target <- function(groups, label_to_target, target_values) {
  # Build a default map: cleaned target value -> target value (identity)
  tl  <- setNames(target_values, clean_label(target_values))
  map <- if (is.null(label_to_target)) tl
  else {
    m <- label_to_target; names(m) <- clean_label(names(m)); m
  }
  # Check every group label has a mapping
  mg <- setdiff(groups, names(map))
  if (length(mg) > 0) stop(
    "Labels in sim_file have no target mapping: ", paste(mg, collapse = ", "), "\n",
    "  These labels came from label1/label2 columns of sim_file (after clean_label).\n",
    "  Targets available in scores_file: ", paste(clean_label(target_values), collapse = ", "), "\n\n",
    "  FIX A: Set LABEL_TO_TARGET to map each label to its target, e.g.:\n",
    "    LABEL_TO_TARGET <- c(\n",
    "      ", paste(sprintf('"%s" = "YOUR_TARGET"', mg), collapse = ",\n      "), "\n",
    "    )\n\n",
    "  FIX B: If labels in sim_file are HPO terms (e.g. HP:0000001), your\n",
    "    sim_file may have the wrong column order. See INPUT FILE FORMATS\n",
    "    at the top of this script.")
  # Check every mapped target exists in scores_file
  mt <- setdiff(unique(unname(map[groups])), target_values)
  if (length(mt) > 0) stop(
    "Mapped targets not found in scores_file: ", paste(mt, collapse = ", "), "\n",
    "  Check spelling and that TARGET_TYPE is correct.")
  map
}

# ---------- New: patient x target matrices (contrastive PhenoSS distance) -
build_patient_target_matrices <- function(phenoss, ids, target_col) {
  all_targets <- sort(unique(as.character(phenoss[[target_col]])))
  score_mat <- matrix(NA_real_, nrow = length(ids), ncol = length(all_targets),
                      dimnames = list(ids, all_targets))
  rank_mat  <- matrix(NA_real_, nrow = length(ids), ncol = length(all_targets),
                      dimnames = list(ids, all_targets))
  for (i in seq_len(nrow(phenoss))) {
    pid <- phenoss$patient_id[i]
    tgt <- as.character(phenoss[[target_col]][i])
    if (pid %in% ids && tgt %in% all_targets) {
      score_mat[pid, tgt] <- phenoss$score[i]
      rank_mat[pid, tgt]  <- phenoss$rank_target[i]
    }
  }
  # Impute: missing ranks -> worst rank for that patient
  for (p in ids) {
    rr <- rank_mat[p, ]
    if (any(is.na(rr))) {
      worst <- if (all(is.na(rr))) length(all_targets) else max(rr, na.rm = TRUE)
      rank_mat[p, is.na(rr)] <- worst
    }
  }
  # Impute: missing scores -> column mean
  for (t in all_targets) {
    sc <- score_mat[, t]
    if (any(is.na(sc))) {
      mn <- if (all(is.na(sc))) 0 else mean(sc, na.rm = TRUE)
      score_mat[is.na(sc), t] <- mn
    }
  }
  # rank_frac per patient: 1 = top rank, 0 = worst rank
  n_cand        <- length(all_targets)
  rank_frac_mat <- if (n_cand == 1)
    matrix(1, nrow = length(ids), ncol = 1, dimnames = list(ids, all_targets))
  else 1 - (rank_mat - 1) / (n_cand - 1)
  list(score = score_mat, rank = rank_mat, rank_frac = rank_frac_mat,
       targets = all_targets)
}

# Contrastive PhenoSS distance:
#   w_p(t)   = -score[p, t] * rank_frac[p, t]    (higher = stronger match)
#   D[i, j] = -(w_i(cj) + w_j(ci))               (lower = patients endorse each other's target)
# Handles vague ranks (e.g. rank 10-50): score magnitude carries the signal
# when rank alone is uninformative.
build_contrastive_dist <- function(ptm, ids, patient_causal) {
  n     <- length(ids)
  w_mat <- (-ptm$score) * ptm$rank_frac   # weight matrix
  D     <- matrix(0, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      ci <- patient_causal[ids[i]]
      cj <- patient_causal[ids[j]]
      if (is.na(ci) || is.na(cj) ||
          !(ci %in% ptm$targets) || !(cj %in% ptm$targets)) {
        D[i, j] <- NA_real_; next
      }
      # i endorses j's target; j endorses i's target
      D[i, j] <- -(w_mat[ids[i], cj] + w_mat[ids[j], ci])
    }
  }
  if (any(is.na(D))) D[is.na(D)] <- max(D, na.rm = TRUE)
  diag(D) <- 0
  D <- (D + t(D)) / 2
  out <- safe_rescale01(D)
  rownames(out) <- ids; colnames(out) <- ids
  out
}

# ---------- Embedding ----------------------------------------------------
embed_2d <- function(dist_mat, method, ids,
                     tsne_perplexity = "auto", tsne_seed = 42,
                     umap_n_neighbors = "auto", umap_min_dist = 0.1,
                     umap_seed = 42) {
  n      <- nrow(dist_mat)
  method <- toupper(method)
  if (method == "MDS") {
    pts <- cmdscale(as.dist(dist_mat), k = 2, eig = TRUE)$points
  } else if (method == "TSNE") {
    perp <- if (identical(tsne_perplexity, "auto"))
      max(1, min(floor(sqrt(n)), floor((n - 1) / 3)))
    else max(1, min(as.numeric(tsne_perplexity), floor((n - 1) / 3)))
    message(sprintf("    t-SNE perplexity = %d  (n=%d)", perp, n))
    set.seed(tsne_seed)
    pts <- Rtsne(as.dist(dist_mat), dims = 2, perplexity = perp,
                 is_distance = TRUE, check_duplicates = FALSE, verbose = FALSE)$Y
  } else if (method == "UMAP") {
    n_nb <- if (identical(umap_n_neighbors, "auto"))
      max(5, min(15, floor(n / 10)))
    else min(as.integer(umap_n_neighbors), n - 1)
    n_nb <- min(n_nb, n - 1)
    message(sprintf("    UMAP n_neighbors = %d  (n=%d)", n_nb, n))
    dist_m  <- as.matrix(dist_mat)
    nn_idx  <- t(apply(dist_m, 1, function(row) order(row)[2:(n_nb + 1)]))
    nn_dist <- t(apply(dist_m, 1, function(row) sort(row)[2:(n_nb + 1)]))
    set.seed(umap_seed)
    pts <- uwot::umap(X = NULL,
                      nn_method    = list(idx = nn_idx, dist = nn_dist),
                      n_components = 2, min_dist = umap_min_dist, verbose = FALSE)
  } else stop("Unknown coord method: ", method, ". Use 'MDS', 'TSNE', or 'UMAP'.")
  rownames(pts) <- ids
  pts
}

# ---------- Alpha label for figure annotations ---------------------------
pub_alpha_label <- function(av) {
  if (av == "resnik_only") return("Resnik only (α=1.0)")
  a <- as.numeric(av)
  if (a == 0.0) return("PhenoSS only (α=0.0)")
  if (a == 0.5) return("Combined ★ (α=0.5)")
  sprintf("α=%.1f", a)
}

# ---------- Static ggplot2 subplot builder -------------------------------
make_static_subplot <- function(pts_2d, groups, dot_sizes, present_groups,
                                color_map, shape_map,
                                overlay_mode, vor,
                                R2_val, p_val, alpha_label,
                                hdb_labels = NULL, km_labels = NULL,
                                hdb_acc = NULL, km_acc = NULL,
                                hdb_sil = NULL, km_sil = NULL,
                                has_true_labels    = TRUE,
                                n_clusters_setting = "true",
                                outlier_iqr        = 3.0,
                                show_legend        = FALSE) {
  stopifnot(length(overlay_mode) == 1,
            overlay_mode %in% c("voronoi","density","hdbscan","kmeans"))

  df   <- data.frame(x = pts_2d[,1], y = pts_2d[,2],
                     group = factor(groups, levels = present_groups),
                     size  = dot_sizes, stringsAsFactors = FALSE)
  keep <- trim_outliers_2d(pts_2d, iqr_thresh = outlier_iqr)
  df   <- df[keep, ]
  group_counts <- table(df$group)
  df$group     <- factor(df$group, levels = names(sort(group_counts, decreasing = TRUE)))
  df           <- df[order(df$group), ]

  xpad <- diff(range(df$x)) * 0.04; ypad <- diff(range(df$y)) * 0.04
  xlim_plot <- range(df$x) + c(-xpad, xpad)
  ylim_plot <- range(df$y) + c(-ypad, ypad)
  leg_pos <- if (show_legend) c(0.98, 0.98) else "none"
  p <- ggplot(df, aes(x = x, y = y, color = group, shape = group, size = size))

  if (overlay_mode == "voronoi") {
    xi <- which(vor$x_seq >= xlim_plot[1] & vor$x_seq <= xlim_plot[2])
    yi <- which(vor$y_seq >= ylim_plot[1] & vor$y_seq <= ylim_plot[2])
    if (length(xi) > 1 && length(yi) > 1) {
      z_sub  <- vor$z_mat[xi, yi, drop = FALSE]
      vor_df <- expand.grid(x = vor$x_seq[xi], y = vor$y_seq[yi])
      vor_df$group <- factor(present_groups[as.vector(z_sub)], levels = present_groups)
      p <- p +
        geom_tile(data = vor_df, aes(x = x, y = y, fill = group),
                  inherit.aes = FALSE, alpha = 0.18, show.legend = FALSE) +
        scale_fill_manual(values = color_map, drop = FALSE)
    }
    cent_df <- as.data.frame(vor$centroids)
    p <- p + geom_point(data = cent_df, aes(x = cx, y = cy),
                        inherit.aes = FALSE, shape = 8, size = 3.5,
                        color = "black", stroke = 1.0)

  } else if (overlay_mode == "density") {
    p <- p + stat_density_2d(aes(color = group), contour = TRUE,
                             linewidth = 0.55, alpha = 0.75, show.legend = FALSE)

  } else if (overlay_mode %in% c("hdbscan", "kmeans")) {
    ul <- if (overlay_mode == "hdbscan") hdb_labels else km_labels
    if (!is.null(ul) && !identical(n_clusters_setting, "true")) {
      cl_df <- data.frame(x = pts_2d[keep, 1], y = pts_2d[keep, 2],
                          cluster = factor(ul[keep]))
      for (cl in levels(cl_df$cluster)) {
        sub <- cl_df[cl_df$cluster == cl, ]
        if (nrow(sub) < 3) next
        hull_idx <- chull(sub$x, sub$y)
        hull_df  <- sub[c(hull_idx, hull_idx[1]), ]
        cl_col   <- make_cluster_color_map(max(as.integer(levels(cl_df$cluster))))[cl]
        p <- p + geom_polygon(data = hull_df, aes(x = x, y = y),
                              inherit.aes = FALSE, fill = cl_col, alpha = 0.12, color = NA)
      }
    }
  }

  p <- p +
    geom_point(alpha = 0.88, stroke = 0.3) +
    scale_color_manual(values = color_map, name = "Group", drop = FALSE) +
    scale_shape_manual(values = shape_map, name = "Group", drop = FALSE) +
    scale_size_identity() +
    coord_cartesian(xlim = xlim_plot, ylim = ylim_plot, expand = FALSE) +
    labs(x = "Dim 1", y = "Dim 2") +
    theme_classic(base_size = 11) +
    theme(
      plot.title       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line        = element_blank(),
      legend.position  = leg_pos,
      legend.justification = c("right","top"),
      legend.background = element_rect(fill = "white", color = "grey60", linewidth = 0.3),
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 10),
      legend.key.size  = unit(0.55, "cm")
    )

  acc_val  <- switch(overlay_mode,
                     voronoi = vor$region_acc, hdbscan = hdb_acc,
                     kmeans = km_acc, density = NULL)
  sil_val  <- switch(overlay_mode, hdbscan = hdb_sil, kmeans = km_sil, NULL)
  acc_type <- if (overlay_mode == "voronoi") "Region acc" else "Cluster acc"

  stat_lines <- c(
    alpha_label,
    sprintf("PERMANOVA R²=%.2f  p=%s", R2_val, signif(p_val, 3)),
    if (!is.null(acc_val) && !is.na(acc_val)) sprintf("%s = %.4f", acc_type, acc_val),
    if (!is.null(sil_val) && !is.na(sil_val)) sprintf("Silhouette = %.4f", sil_val),
    "● size = PhenoSS rank"
  )
  p + annotate("label",
               x = xlim_plot[2], y = ylim_plot[1],
               label = paste(stat_lines, collapse = "\n"),
               hjust = 1.0, vjust = 0.0, size = 2.5, lineheight = 1.25,
               label.padding = unit(0.45, "lines"), label.r = unit(0.12, "lines"),
               fill = "white", color = "black", fontface = "plain", alpha = 0.93)
}

# ---------- Interactive Plotly builder (used by commented-out STEP 10) ---
make_plot <- function(hover_data, title_str, vor,
                      present_groups, group_colors, color_map, symbol_map,
                      R2_val, p_val, overlay_mode, alpha_label, coord_method,
                      has_true_labels = TRUE, unsup_labels = NULL,
                      cluster_color_map = NULL, accuracy_val = NULL,
                      silhouette_val = NULL, n_clusters_setting = "true") {
  fig    <- plotly::plot_ly()
  shapes <- list()
  if (overlay_mode == "voronoi") {
    rects  <- make_voronoi_rects(vor$x_seq, vor$y_seq, vor$z_mat,
                                 present_groups, group_colors, alpha = 0.20, ds = 6)
    bnd    <- get_boundary_shapes(vor$x_seq, vor$y_seq, vor$z_mat)
    shapes <- c(shapes, rects, bnd)
    fig    <- plotly::add_markers(fig, x = vor$centroids$cx, y = vor$centroids$cy,
               marker = list(symbol = "star", size = 16, color = "black",
                             line = list(color = "white", width = 1.5)),
               name = "Centroids", hoverinfo = "skip", showlegend = TRUE)
  } else if (overlay_mode == "density") {
    for (g in present_groups) {
      sub <- hover_data[hover_data$group == g, ]
      if (nrow(sub) < 3) next
      fig <- plotly::add_trace(fig, x = sub$dim1, y = sub$dim2,
               type = "histogram2dcontour", ncontours = 5, showscale = FALSE,
               contours = list(coloring = "lines"),
               line = list(color = color_map[g], width = 1.5),
               hoverinfo = "skip", showlegend = FALSE)
    }
  } else if (overlay_mode %in% c("hdbscan","kmeans") &&
             !is.null(unsup_labels) && !identical(n_clusters_setting, "true")) {
    for (cl in sort(unique(unsup_labels))) {
      idx <- which(unsup_labels == cl)
      if (length(idx) < 3) next
      pts <- cbind(hover_data$dim1[idx], hover_data$dim2[idx])
      hi  <- chull(pts); hi <- c(hi, hi[1])
      col  <- cluster_color_map[as.character(cl)]; rgbv <- col2rgb(col)
      fillc <- sprintf("rgba(%d,%d,%d,0.12)", rgbv[1], rgbv[2], rgbv[3])
      path  <- paste0("M ", pts[hi[1],1], ",", pts[hi[1],2], " ",
                      paste(sprintf("L %s,%s", pts[hi[-1],1], pts[hi[-1],2]), collapse=" "), " Z")
      shapes[[length(shapes)+1]] <- list(type="path", layer="below", path=path,
        fillcolor=fillc, line=list(color=col, width=1))
    }
  }
  for (g in present_groups) {
    sub <- hover_data[hover_data$group == g, ]
    if (nrow(sub) == 0) next
    fig <- plotly::add_markers(fig, x=sub$dim1, y=sub$dim2,
             text=sub$hover_text, hoverinfo="text", name=g,
             marker=list(symbol=symbol_map[g], color=color_map[g],
                         size=sub$dot_size, line=list(color="black", width=0.5)))
  }
  acc_txt <- if (!is.null(accuracy_val) && !is.na(accuracy_val)) {
    lbl <- if (overlay_mode == "voronoi") "Region acc" else "Cluster acc"
    sprintf("<br>%s = %.4f", lbl, accuracy_val)
  } else ""
  sil_txt <- if (!is.null(silhouette_val) && !is.na(silhouette_val))
    sprintf("<br>Silhouette = %.4f", silhouette_val) else ""
  ann_txt <- sprintf("<b>%s</b><br>PERMANOVA R²=%.2f  p=%s%s%s",
                     alpha_label, R2_val, signif(p_val,3), acc_txt, sil_txt)
  plotly::layout(fig,
    title  = list(text = title_str, x = 0.5),
    xaxis  = list(title = "Dim 1", zeroline = FALSE),
    yaxis  = list(title = "Dim 2", zeroline = FALSE),
    shapes = shapes,
    annotations = list(list(x=1, y=0, xref="paper", yref="paper",
      xanchor="right", yanchor="bottom", text=ann_txt, showarrow=FALSE,
      bgcolor="rgba(255,255,255,0.9)", bordercolor="grey", borderwidth=1,
      borderpad=6, font=list(size=12))),
    legend = list(x=1.02, y=1, bgcolor="rgba(255,255,255,0.8)"))
}

# =============================================================================
# STEP 1 — SETUP
# =============================================================================
TARGET_COL <- resolve_target_column(TARGET_TYPE)
message("TARGET_TYPE=", TARGET_TYPE, "  column=", TARGET_COL)
if (!dir.exists(VIZ_DIR)) { dir.create(VIZ_DIR, recursive = TRUE); message("Created: ", VIZ_DIR) }

# =============================================================================
# STEP 2 — RESNIK SIMILARITY FILE
#
# Expected format (tab-separated, NO header):
#   col 1: patient_id_1
#   col 2: patient_id_2
#   col 3: similarity (numeric > 0)
#   col 4: label for patient 1 (gene/disease group name)
#   col 5: label for patient 2 (gene/disease group name)
#
# If your file has a header row, change header=FALSE to header=TRUE below
# and remove the colnames() assignment.
# =============================================================================
x <- read.table(sim_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(x) <- c("id1", "id2", "sim", "label1", "label2")
x <- x[!is.na(x$sim), ]
x$id1 <- trimws(x$id1); x$id2 <- trimws(x$id2)

# --- Preflight: catch wrong column order early ---------------------------
# This fires if sim_file has HPO terms in column 3 instead of a number,
# which happens when phenoss_similarity.tsv (col order: pat1 pat2 hpo1 hpo2 sim)
# is used directly without reordering columns first.
sim_numeric <- suppressWarnings(as.numeric(x$sim))
if (all(is.na(sim_numeric))) stop(
  "Column 3 of sim_file is not numeric. Got: '", x$sim[1], "'\n\n",
  "  This usually means you are using phenoss_similarity.tsv directly,\n",
  "  which has column order: pat1  pat2  hpo1  hpo2  similarity\n",
  "  but clustering.R expects:  id1   id2   sim   label1  label2\n\n",
  "  You need to reorder columns AND replace hpo1/hpo2 with group labels.\n",
  "  See INPUT FILE FORMATS at the top of this script for a Python snippet.")
x$sim <- sim_numeric

# Capture original (pretty) label casing before clean_label() lowercases them.
# Used later by prettify_label() in figure titles.
raw_pretty <- unique(c(trimws(as.character(x$label1)), trimws(as.character(x$label2))))
PRETTY_LABEL_MAP <- setNames(raw_pretty, clean_label(raw_pretty))

x$label1 <- clean_label(x$label1); x$label2 <- clean_label(x$label2)
x$dist   <- pmin(1 / x$sim, 1)
stopifnot(all(x$dist > 0), all(x$dist <= 1), all(is.finite(x$dist)))

# =============================================================================
# STEP 3 — PHENOSS SCORES FILE
#
# This is the output from run_phenoss.py.
# Required columns: patient_id, score, and disease_mondo (or gene).
# =============================================================================
scores_raw    <- read.table(scores_file, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, quote = "")
phenoss       <- prepare_scores(scores_raw, TARGET_COL, TARGET_TYPE)
target_values <- sort(unique(as.character(phenoss[[TARGET_COL]])))
message("Targets found in scores_file: ", paste(target_values, collapse = ", "))

# =============================================================================
# STEP 4 — PATIENT IDs + GROUP METADATA
# =============================================================================
ids        <- sort(unique(c(x$id1, x$id2)))
n          <- length(ids)
raw_groups <- sort(unique(c(x$label1, x$label2)))

has_true_labels_global <- length(raw_groups) > 0 &&
  !all(raw_groups %in% c("", "na", "unknown"))
n_true_labels <- if (has_true_labels_global) length(raw_groups) else 0L

# Build label -> target mapping.
# If LABEL_TO_TARGET = NULL: expects labels to already match target IDs.
# Otherwise: maps each sim_file label to its target ID in scores_file.
label_to_target_map <- resolve_label_to_target(raw_groups, LABEL_TO_TARGET, target_values)

if (is.null(GROUP_ORDER)) GROUP_ORDER <- sort(raw_groups)
meta           <- build_meta(x, ids, group_order = GROUP_ORDER)
present_groups <- levels(droplevels(meta$group))
true_int       <- as.integer(meta$group)

message(sprintf("Patients: %d  |  Groups (%d): %s",
                n, n_true_labels, paste(present_groups, collapse = ", ")))

k_resolved <- resolve_k(N_CLUSTERS, n_true_labels)
message(sprintf("N_CLUSTERS='%s' -> k=%s",
                N_CLUSTERS, if (is.na(k_resolved)) "auto" else k_resolved))

# =============================================================================
# STEP 5 — COLORS + SHAPES
# =============================================================================
color_map    <- make_color_map(present_groups)
shape_map    <- make_shape_map(present_groups)
symbol_map   <- make_symbol_map(present_groups)
group_colors <- color_map[present_groups]

# =============================================================================
# STEP 6 — RESNIK DISTANCE MATRIX  (from sim_file similarity scores)
# =============================================================================
dist_mat_pheno <- matrix(1, n, n, dimnames = list(ids, ids))
diag(dist_mat_pheno) <- 0
for (i in seq_len(nrow(x))) {
  dist_mat_pheno[x$id1[i], x$id2[i]] <- x$dist[i]
  dist_mat_pheno[x$id2[i], x$id1[i]] <- x$dist[i]
}

# =============================================================================
# STEP 7 — PHENOSS TARGET-SCORE DISTANCE MATRIX
#
# Uses a cross-patient contrastive formulation across ALL causal targets
# in the cohort (not just each patient's own target scalar).
#
# For patients i and j with causal targets ci and cj:
#   w_p(t)   = -score[p, t] * rank_frac[p, t]    (higher = stronger match)
#   D[i, j] = -(w_i(cj) + w_j(ci))
#
# This means: two patients are close when each strongly endorses the other's
# causal target — capturing cross-disease discrimination even when all patients
# happen to rank their own target at #1.
#
# rank_frac normalizes by total candidates per patient, so patients scored
# against different-sized candidate lists are still comparable.
# =============================================================================
ptm <- build_patient_target_matrices(phenoss, ids, TARGET_COL)
message(sprintf("Score matrix: %d patients x %d targets",
                nrow(ptm$score), ncol(ptm$score)))

patient_causal <- unname(label_to_target_map[as.character(meta$group)])
names(patient_causal) <- ids

dist_mat_target <- build_contrastive_dist(ptm, ids, patient_causal)
message("Contrastive PhenoSS distance matrix built.")

# Keep per-patient stats for dot sizes and hover text
chosen_stats <- data.frame(
  patient_id   = ids,
  target_label = as.character(meta$group),
  target_value = unname(patient_causal[ids]),
  score        = ptm$score[cbind(ids, unname(patient_causal[ids]))],
  rank_target  = ptm$rank[ cbind(ids, unname(patient_causal[ids]))],
  rank_frac    = ptm$rank_frac[cbind(ids, unname(patient_causal[ids]))],
  stringsAsFactors = FALSE
)

# =============================================================================
# STEP 8 — DOT SIZES  (larger dot = better PhenoSS rank for own causal target)
# =============================================================================
proband_rank <- setNames(chosen_stats$rank_target, chosen_stats$patient_id)[ids]
dot_sizes_gg <- make_dot_sizes(proband_rank, min_size = 1.5, max_size = 5)

# =============================================================================
# STEP 9 — PERMANOVA + FUSED DISTANCE MATRICES
#
# For each alpha:
#   D_combined = alpha * D_resnik + (1 - alpha) * D_phenoss
#   alpha = 1.0  -> Resnik only
#   alpha = 0.0  -> PhenoSS only (contrastive score x rank_frac)
#   alpha = 0.5  -> equal blend
#
# PERMANOVA R² tells you how much variance in distances is explained by group.
# Higher R² and lower p = better separation by label.
# =============================================================================
message("\n=== Pre-computing PERMANOVA stats ===")
perm_stats <- list()

message("  resnik_only ...")
pr <- adonis2(as.dist(dist_mat_pheno) ~ group, data = meta, permutations = 999)
perm_stats[["resnik_only"]] <- list(R2  = round(pr$R2[1], 2),
                                    p   = pr$`Pr(>F)`[1],
                                    mat = dist_mat_pheno)

for (av in ALPHA_VALUES[ALPHA_VALUES != "resnik_only"]) {
  a   <- as.numeric(av)
  key <- paste0("alpha_", a)
  # Combined distance: alpha=1 -> Resnik only; alpha=0 -> PhenoSS only
  mat <- a * dist_mat_pheno + (1 - a) * dist_mat_target
  diag(mat) <- 0; mat <- (mat + t(mat)) / 2
  message(sprintf("  alpha=%.2f ...", a))
  pa  <- adonis2(as.dist(mat) ~ group, data = meta, permutations = 999)
  perm_stats[[key]] <- list(R2  = round(pa$R2[1], 2),
                             p   = pa$`Pr(>F)`[1],
                             mat = mat)
}
message("  Done. R² summary:")
for (k in names(perm_stats))
  message(sprintf("    %-20s  R²=%.2f  p=%.4f", k, perm_stats[[k]]$R2, perm_stats[[k]]$p))

# =============================================================================
# STEP 10 — INTERACTIVE PLOTLY LOOP  (commented out; preserved for future use)
#
# To enable:
#   1. Uncomment library(plotly) and library(htmlwidgets) at the top.
#   2. Uncomment this entire block.
# =============================================================================
# total_plots  <- length(ALPHA_VALUES) * length(COORD_METHODS) * length(OVERLAY_MODES)
# plot_counter <- 0
# message(sprintf("\n=== Generating %d interactive HTML plots ===\n", total_plots))
#
# for (av in ALPHA_VALUES) {
#   key          <- if (av == "resnik_only") "resnik_only" else paste0("alpha_", as.numeric(av))
#   alpha_label  <- if (av == "resnik_only") "resnik_only" else sprintf("alpha%.2f", as.numeric(av))
#   R2_val       <- perm_stats[[key]]$R2
#   p_val        <- perm_stats[[key]]$p
#   dist_mat_use <- perm_stats[[key]]$mat
#
#   for (cm in COORD_METHODS) {
#     message(sprintf("  Embedding: %s | %s", alpha_label, cm))
#     pts_2d <- tryCatch(
#       embed_2d(dist_mat_use, cm, ids,
#                tsne_perplexity=TSNE_PERPLEXITY, tsne_seed=TSNE_SEED,
#                umap_n_neighbors=UMAP_N_NEIGHBORS, umap_min_dist=UMAP_MIN_DIST,
#                umap_seed=UMAP_SEED),
#       error=function(e){message("  ERROR: ",e$message); NULL})
#     if (is.null(pts_2d)) next
#
#     vor              <- compute_voronoi_objects(pts_2d, meta$group, present_groups,
#                                                GRID_RES, PAD_FRAC, GROUP_BIAS)
#     dot_sizes_plotly <- make_dot_sizes(proband_rank, min_size=7, max_size=22)
#     tln              <- if (tolower(TARGET_TYPE)=="gene") "Target gene" else "Target disease"
#     hover_df <- data.frame(
#       id=ids, group=as.character(meta$group),
#       target_value=chosen_stats$target_value[match(ids,chosen_stats$patient_id)],
#       rank=as.numeric(proband_rank),
#       score=chosen_stats$score[match(ids,chosen_stats$patient_id)],
#       dim1=pts_2d[,1], dim2=pts_2d[,2], dot_size=dot_sizes_plotly,
#       stringsAsFactors=FALSE)
#     hover_df$hover_text <- paste0(
#       "<b>",hover_df$id,"</b><br>Label: ",hover_df$group,"<br>",
#       tln,": ",hover_df$target_value,"<br>PhenoSS rank: ",hover_df$rank,
#       "<br>Score: ",round(hover_df$score,6))
#
#     hdb_labels <- tryCatch(run_hdbscan(pts_2d,k_resolved),error=function(e)NULL)
#     km_labels  <- tryCatch(run_kmeans(pts_2d,k_resolved), error=function(e)NULL)
#     hdb_acc    <- if(n_true_labels>0 && !is.null(hdb_labels)) compute_cluster_accuracy(true_int,hdb_labels) else NULL
#     km_acc     <- if(n_true_labels>0 && !is.null(km_labels))  compute_cluster_accuracy(true_int,km_labels)  else NULL
#     hdb_sil    <- if(!is.null(hdb_labels)) compute_silhouette_score(pts_2d,hdb_labels) else NA
#     km_sil     <- if(!is.null(km_labels))  compute_silhouette_score(pts_2d,km_labels)  else NA
#     k_n        <- if(is.na(k_resolved)) max(2L,length(unique(c(hdb_labels,km_labels))),na.rm=TRUE) else k_resolved
#     cl_colors  <- make_cluster_color_map(k_n)
#
#     for (om in OVERLAY_MODES) {
#       plot_counter <- plot_counter + 1
#       fname <- sprintf("%s_%s_%s_%s.html", CASE_NAME, cm, alpha_label, om)
#       fpath <- file.path(VIZ_DIR, fname)
#       message(sprintf("  [%d/%d] %s", plot_counter, total_plots, fname))
#       ul  <- if(om=="hdbscan") hdb_labels else if(om=="kmeans") km_labels else NULL
#       acc <- if(om=="hdbscan") hdb_acc    else if(om=="kmeans") km_acc    else NULL
#       sil <- if(om=="hdbscan") hdb_sil    else if(om=="kmeans") km_sil    else NULL
#       fig <- tryCatch(
#         make_plot(hover_data=hover_df, title_str=sprintf("%s | %s | %s | %s", CASE_NAME, cm, alpha_label, om),
#                   vor=vor, present_groups=present_groups, group_colors=group_colors,
#                   color_map=color_map, symbol_map=symbol_map,
#                   R2_val=R2_val, p_val=p_val, overlay_mode=om,
#                   alpha_label=alpha_label, coord_method=cm,
#                   has_true_labels=(n_true_labels>0), unsup_labels=ul,
#                   cluster_color_map=cl_colors, accuracy_val=acc,
#                   silhouette_val=sil, n_clusters_setting=N_CLUSTERS),
#         error=function(e){message("  ERROR: ",e$message); NULL})
#       if (!is.null(fig)) {
#         tmp <- tempfile(fileext=".html")
#         htmlwidgets::saveWidget(fig, tmp, selfcontained=TRUE)
#         file.copy(tmp, fpath, overwrite=TRUE); file.remove(tmp)
#         message("    Saved -> ", fpath)
#       }
#     }
#   }
# }

# =============================================================================
# STEP 11 — STATIC PUBLICATION FIGURES
#
# Produces one figure per (coord method, overlay combination).
# Each figure has len(PUB_ALPHAS) panels side-by-side.
# Output: PDF + TIFF + PNG per figure.
#
# Filenames: {CASE_NAME}_publication_{coord}_{overlay}.{pdf|tiff|png}
# =============================================================================
message("\n=== Generating static publication figures ===")
message(sprintf("  %d coord methods x %d overlays = %d figures",
                length(COORD_METHODS), length(PUB_OVERLAYS),
                length(COORD_METHODS) * length(PUB_OVERLAYS)))

for (cm in COORD_METHODS) {
  message(sprintf("\n  Coord method: %s", cm))

  # Compute embeddings once per (coord, alpha); share across overlays
  embed_cache <- vector("list", length(PUB_ALPHAS))
  for (i in seq_along(PUB_ALPHAS)) {
    av  <- PUB_ALPHAS[i]
    key <- if (av == "resnik_only") "resnik_only" else paste0("alpha_", as.numeric(av))
    if (is.null(perm_stats[[key]])) {
      warning("No perm_stats for '", key, "' — skipping. ",
              "Add it to ALPHA_VALUES to compute it."); next
    }
    dist_mat_use <- perm_stats[[key]]$mat
    message(sprintf("    Embedding: %s", pub_alpha_label(av)))
    pts_2d <- tryCatch(
      embed_2d(dist_mat_use, cm, ids,
               tsne_perplexity  = TSNE_PERPLEXITY,  tsne_seed = TSNE_SEED,
               umap_n_neighbors = UMAP_N_NEIGHBORS,  umap_min_dist = UMAP_MIN_DIST,
               umap_seed = UMAP_SEED),
      error = function(e) { message("    ERROR embedding: ", e$message); NULL }
    )
    if (is.null(pts_2d)) next
    vor        <- compute_voronoi_objects(pts_2d, meta$group, present_groups,
                                          GRID_RES, PAD_FRAC, GROUP_BIAS)
    hdb_labels <- tryCatch(run_hdbscan(pts_2d, k_resolved), error = function(e) NULL)
    km_labels  <- tryCatch(run_kmeans(pts_2d, k_resolved),  error = function(e) NULL)
    has_true   <- n_true_labels > 0
    hdb_acc <- if (has_true && !is.null(hdb_labels)) compute_cluster_accuracy(true_int, hdb_labels) else NULL
    km_acc  <- if (has_true && !is.null(km_labels))  compute_cluster_accuracy(true_int, km_labels)  else NULL
    hdb_sil <- if (!is.null(hdb_labels)) compute_silhouette_score(pts_2d, hdb_labels) else NA
    km_sil  <- if (!is.null(km_labels))  compute_silhouette_score(pts_2d, km_labels)  else NA
    embed_cache[[i]] <- list(
      pts_2d = pts_2d, vor = vor,
      hdb_labels = hdb_labels, km_labels = km_labels,
      hdb_acc = hdb_acc, km_acc = km_acc,
      hdb_sil = hdb_sil, km_sil = km_sil,
      R2_val = perm_stats[[key]]$R2, p_val = perm_stats[[key]]$p,
      alabel = pub_alpha_label(av), has_true = has_true
    )
  }

  # One figure per overlay, reusing the cached embeddings
  for (om in PUB_OVERLAYS) {
    message(sprintf("    Overlay: %s", om))
    subplots <- vector("list", length(PUB_ALPHAS))

    for (i in seq_along(PUB_ALPHAS)) {
      ec <- embed_cache[[i]]
      if (is.null(ec)) next
      subplots[[i]] <- make_static_subplot(
        pts_2d             = ec$pts_2d,
        groups             = as.character(meta$group),
        dot_sizes          = dot_sizes_gg,
        present_groups     = present_groups,
        color_map          = color_map,
        shape_map          = shape_map,
        overlay_mode       = om,
        vor                = ec$vor,
        R2_val             = ec$R2_val,
        p_val              = ec$p_val,
        alpha_label        = ec$alabel,
        hdb_labels         = ec$hdb_labels,
        km_labels          = ec$km_labels,
        hdb_acc            = ec$hdb_acc,
        km_acc             = ec$km_acc,
        hdb_sil            = ec$hdb_sil,
        km_sil             = ec$km_sil,
        has_true_labels    = ec$has_true,
        n_clusters_setting = N_CLUSTERS,
        outlier_iqr        = OUTLIER_IQR_THRESHOLD,
        show_legend        = (i == length(PUB_ALPHAS))
      )
    }

    subplots <- Filter(Negate(is.null), subplots)
    if (length(subplots) == 0) { message("    No panels — skipping."); next }

    # Dynamic title from the actual group labels in sim_file
    cohort_str <- format_label_list(vapply(present_groups, prettify_label, character(1)))
    combined <- wrap_plots(subplots, nrow = 1) +
      plot_annotation(
        title = sprintf("Patient Clustering of %s Syndromes  |  %s  |  %s",
                        cohort_str, cm, om),
        theme = theme(
          plot.title      = element_text(size = 15, face = "bold",
                                         hjust = 0.5, margin = margin(b = 10)),
          plot.background = element_rect(fill = "white", color = NA)
        )
      )

    base      <- sprintf("%s_publication_%s_%s", CASE_NAME, cm, om)
    pdf_path  <- file.path(VIZ_DIR, paste0(base, ".pdf"))
    tiff_path <- file.path(VIZ_DIR, paste0(base, ".tiff"))
    png_path  <- file.path(VIZ_DIR, paste0(base, ".png"))
    ggsave(pdf_path,  combined, width = PUB_WIDTH, height = PUB_HEIGHT,
           units = "in", device = cairo_pdf)
    ggsave(tiff_path, combined, width = PUB_WIDTH, height = PUB_HEIGHT,
           units = "in", dpi = PUB_DPI, device = "tiff", compression = "lzw")
    ggsave(png_path,  combined, width = PUB_WIDTH, height = PUB_HEIGHT,
           units = "in", dpi = PUB_DPI)
    message(sprintf("      Saved: %s.{pdf,tiff,png}", base))
  }
}

message("\n=== Done. Figures saved to: ", VIZ_DIR, " ===")
