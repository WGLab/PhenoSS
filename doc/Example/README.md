# PhenoSS Tutorial
## Overview
This tutorial demonstrates how to use PhenoSS to perform patient clustering. The files are contained the folder Example.

## Patient Clustering
### Sample HPO data
The file hpo_list contains the synthetic data for three randomly generated patients labeled 0_10, 1_10, 2_10. 
```
0_10    HP_0004370;HP_0000280;HP_0002835;HP_0005274;HP_0000158;HP_0011470;HP_0001417;HP_0001270;HP_0008872;HP_0002015;HP_0000750;HP_0000157;
1_10    HP_0000483;HP_0002307;HP_0001090;HP_0001572;HP_0002342;HP_0011343;HP_0008760;HP_0001061;HP_0001249;HP_0000574;HP_0001417;HP_0002020;HP_0012810;HP_000
0540;HP_0001350;HP_0001270;HP_0002574;HP_0011231;HP_0000750;HP_0002155;HP_0000431;HP_0000718;
2_10    HP_0002194;HP_0001263;HP_0001684;HP_0001417;HP_0001270;HP_0001249;HP_0001670;HP_0001667;HP_0001629;HP_0001639;HP_0002474;HP_0010863;HP_0000750;
```

### Similairty score calculation

Using the following argument, we can calculate the similarity scores between patient 1_10 and each of the patients in the hpo_list. 
The first input argument is the file that contains the patient IDs and the HPO terms. The second input argument is the patient ID we are interested in.
```
python getdiff_one_patient.py hpo_list 1_10
```

You are recommended to submit one job for each patient to perform computing simultaneously.

The outputs of the argument can be found in the file 1_10_sim.
```
1_10	0_10	3.741497748992082
1_10	1_10	8.089540938721928
1_10	2_10	2.451702686307226
```
The following code will summarize the scores for all the patients and form the similarity matrix.
```
perl get_sim_mat.pl hpo_list
```
We can then perform the quality check. filter.pl automatically detects missing values in the similarity matrix and selects te largest subset of the patients such that there is no missing similarity scores for these patients and forms the similarity score matrix.

```
perl filter.pl
```
The similarity score matrix will then be saved into the file `sim_mat_filter`.

### Hierarchical clustering
To better understand the results, we can perform the hierarchical clustering and plot the dendrograms in R.
```
library(dplyr)
x <- read.table("sim_mat_filter")
x$V3 <- 1/x$V3
numpat <- sqrt(dim(x)[1])
pat_mat <- matrix(x$V3, nrow = numpat)
patid <- x$V2[1:numpat]
colnames(pat_mat) <- patid
rownames(pat_mat) <- patid

hc <- hclust(dist(pat_mat), method = "ward.D2")
dend <- as.dendrogram(hc)

#groupCodes <- c(rep("NA10", 63), rep("NA15", 64))
#colorCodes <- c(NA10="blue", NA15="green")
#labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]

plot(dend)
```

