library(dplyr)
x <- read.table("/Users/gracsh/Desktop/EHRclust-master/sim_mat_filter")
x$V3 <- 1/x$V3
numpat <- sqrt(dim(x)[1])
pat_mat <- matrix(x$V3, nrow = numpat)
patid <- x$V2[1:numpat]
colnames(pat_mat) <- patid
rownames(pat_mat) <- patid

hc <- hclust(dist(pat_mat), method = "ward.D2")
dend <- as.dendrogram(hc)

groupCodes <- c(rep("NA10", 1), rep("NA15", 2), rep("NA20", 1))
colorCodes <- c(NA10="blue", NA15="green", NA20="red")
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]

plot(dend)