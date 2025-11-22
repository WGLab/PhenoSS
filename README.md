# PhenoSS: Phenotype semantic similarity-based approach for rare disease prediction and patient clustering 

## Introduction
PhenoSS is an effective algorithm that makes disease prediction and performs patient clustering based on HPO concepts. PhenoSS uses the Gaussian copula technique by modeling the marginal prevalence of each HPO term for each disease and utilizes a multivariate normal distribution to link them together to account for term correlations. We utilized the OARD (open annotations for rare diseases) API for inferring the frequency of HPO terms in a diverse range of rare diseases. PhenoSS can calculate the phenotype similarity between any two patients for finding similar patients or clustering purposes, or between one patient and any candidate diseases for diagnosis support. 

The toolkit is implemented in Python. 

## Installation



## Tutorial
### Patient Clustering
#### Sample HPO data
The file hpo_list contains the synthetic data for three randomly generated patients labeled 0_10, 1_10, 2_10. 
```
0_10    HP_0004370;HP_0000280;HP_0002835;HP_0005274;HP_0000158;HP_0011470;HP_0001417;HP_0001270;HP_0008872;HP_0002015;HP_0000750;HP_0000157;
1_10    HP_0000483;HP_0002307;HP_0001090;HP_0001572;HP_0002342;HP_0011343;HP_0008760;HP_0001061;HP_0001249;HP_0000574;HP_0001417;HP_0002020;HP_0012810;HP_000
0540;HP_0001350;HP_0001270;HP_0002574;HP_0011231;HP_0000750;HP_0002155;HP_0000431;HP_0000718;
2_10    HP_0002194;HP_0001263;HP_0001684;HP_0001417;HP_0001270;HP_0001249;HP_0001670;HP_0001667;HP_0001629;HP_0001639;HP_0002474;HP_0010863;HP_0000750;
```

#### Similairty score calculation

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

#### Hierarchical clustering
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

### Disease prediction
PhenoSS extracts the diseases/phenotype frequencies from the Open Annotations for Rare Diseases (OARD) Database. It takes in HPO terms of a list of patients and outputs the ranks of possible underlying diseases. 

Below is a sample input file:

```
P1	HP_0012759;HP_0000750;HP_0100022;HP_0000707;
P2	HP_0001270;HP_0012758;HP_0002066;HP_0011443;
P3	HP_0012758;HP_0002167;HP_0012638;HP_0000707;
```
To run PhenoSS, use the following command:

```
python phenoSS.py inputFile outputFile
```
The results consist of a list of MONDO diseases and the rankings and will be stored in 'outputFile' specified by the user. 

### Convert MONDO diseases to genes
PhenoSS outputs a list of MONDO diseases and the corresponding rankings. The file 'mondo2gene.txt' maps MONDO diseases to gene symbols. Below is the first four lines of 'mondo2gene.txt':

| class  | class_label | OMIM | Approved Gene Symbol (HGNC) |
| ------------- | ------------- | ------------- | ------------- |
| MONDO:0013138  | BRV2  | OMIM:613106  | BRV2  |
| MONDO:0014068  | cone-rod dystrophy 17  | OMIM:615163  | CORD17  |
| MONDO:0013151  | CACD3  | OMIM:613144  | CACD3  |
| MONDO:0010568  | Aicardi syndrome  | OMIM:304050  | AIC  |


To convert the results into genes, run the following command:

```
python mondo2gene.py inputFile outputFile
```
The file 'inputFile' should contain the output from PhenoSS, and the converted genes will be stored in 'outputFile'.
## Datasets
#### Human Phenotype Ontology (HPO): 
https://hpo.jax.org/
#### Open Annotations for Rare Diseases (OARD):
https://rare.cohd.io/

## License

PhenoSS is distributed under the [MIT License by Wang Genomics Lab](https://wglab.mit-license.org/).
