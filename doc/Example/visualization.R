!install.packages("dendextend")   # run once if not installed
library(dendextend)
library(dplyr)
x <- read.table("/home/nguyenqm/projects/github/PhenoSS/doc/paper_data/sim_mat_deidentified_filter")
x$V3 <- 1/x$V3
numpat <- sqrt(dim(x)[1])
pat_mat <- matrix(x$V3, nrow = numpat)
patid <- x$V2[1:numpat]
colnames(pat_mat) <- patid
rownames(pat_mat) <- patid

icd_code <- sub("_.*", "", patid)
names(icd_code) <- patid
group_map <- c(
  "Q85.01" = "NA20",
  "G11.11" = "NA10",
  "Q87.40" = "NA15"
)
disease_labels <- c(
  NA20 = "Neurofibromatosis type 1",
  NA10 = "Friedreich Ataxia",
  NA15 = "Marfan Syndrome"
)

hc <- hclust(dist(pat_mat), method = "ward.D2")
dend <- as.dendrogram(hc)

#groupCodes <- c(rep("NA10", 1), rep("NA15", 2), rep("NA20", 1))
groupCodes <- group_map[icd_code]
names(groupCodes) <- patid

colorCodes <- c(NA10="blue", NA15="green", NA20="red")
#labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
lab <- labels(dend)

labels_colors(dend) <- colorCodes[groupCodes[lab]]

png("/home/nguyenqm/projects/github/PhenoSS/doc/paper_data/dendrogram_clarity.png", width = 2400, height = 1600, res = 150)

plot(
  dend,
  cex = 0.4,     # smaller label size
  cex.main = 1.5,
  main = "Clustering of Patients with Friedreich Ataxia, Neurofibromatosis Type 1, or Marfan Syndrome"
)

legend(
  "topleft",
  legend = disease_labels[names(colorCodes)],
  fill   = colorCodes,
  border = NA,
  bty    = "n",
  cex    = 1.5
)

dev.off()

