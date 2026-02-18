!install.packages("FNN")
!install.packages("vegan")
library(vegan)
library(FNN)
library(dplyr)

# -----------------------------
# 1. Read similarity table
# -----------------------------
file_path <- "/home/nguyenqm/projects/github/PhenoSS/doc/paper_data/sim_mat_arcus"
x <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)

colnames(x) <- c("id1","id2","sim")

# convert similarity to distance
x$dist <- 1 / x$sim

# -----------------------------
# 2. Build distance matrix
# -----------------------------
ids <- sort(unique(c(x$id1, x$id2)))
n <- length(ids)

mat <- matrix(0, nrow = n, ncol = n)
rownames(mat) <- ids
colnames(mat) <- ids

for(i in seq_len(nrow(x))){
  mat[x$id1[i], x$id2[i]] <- x$dist[i]
}

# make symmetric
mat <- (mat + t(mat))/2

# -----------------------------
# 3. Run MDS
# -----------------------------
mds <- cmdscale(as.dist(mat), k = 2)

# -----------------------------
# 4. Extract disease from ICD code
# -----------------------------
icd <- sub("_.*","", rownames(mds))

group <- ifelse(icd == "G11.11", "Ataxia",
                ifelse(icd == "Q87.40", "Marfan",
                       ifelse(icd == "Q85.01", "Neurofibromatosis T1", "Other")))

colors <- c(
  "Ataxia" = "blue",
  "Marfan" = "green",
  "Neurofibromatosis T1" = "red"
)

# -----------------------------
# PERMANOVA
# -----------------------------
dist_obj <- as.dist(mat)
meta <- data.frame(id = ids, group = group)
perm <- adonis2(dist_obj ~ group, data = meta, permutations = 999)

R2_val <- round(perm$R2[1], 2)
p_val  <- perm$`Pr(>F)`[1]

# -----------------------------
# kNN accuracy (already computed)
# -----------------------------
D <- as.matrix(mat)
y <- factor(group)

# convert distance matrix to neighbor order
# for each i, neighbors are smallest distances excluding itself
pred1 <- sapply(seq_len(n), function(i){
  nn <- order(D[i, ])[2]         # 1-NN
  as.character(y[nn])
})
acc <- mean(pred1 == y)
acc_text <- paste0("1-NN accuracy: ", round(acc, 2))

# -----------------------------
# Plot with metrics
# -----------------------------
# define plotting limits excluding extreme outliers
xlim_use <- quantile(mds[,1], c(0.02, 0.98))
ylim_use <- quantile(mds[,2], c(0.02, 0.98))

png("/home/nguyenqm/projects/github/PhenoSS/doc/paper_data/arcus_mds_plot.png",
    width = 1800, height = 1400, res = 200)

plot(
  mds,
  col = colors[group],
  pch = 1,
  xlab = "Dim 1",
  ylab = "Dim 2",
  xlim = xlim_use,
  ylim = ylim_use,
  main = "MDS of ARCUS Patients with Friedreich Ataxia,\nNeurofibromatosis Type 1, or Marfan Syndrome"
)

legend(
  "bottomright",
  legend = c("Ataxia","Marfan","Neurofibromatosis T1"),
  col = colors,
  pch = 15,
  bty = "n"
)

text(
  x = xlim_use[1],
  y = ylim_use[1],
  adj = c(0,0),
  labels = paste0(
    "PERMANOVA: R² = ", R2_val,
    ", p = ", p_val,
    "\n", acc_text
  ),
  cex = 1.1
)

dev.off()

