# PhenoSS: Phenotype semantic similarity-based approach for rare disease prediction and patient clustering 

# Introduction
PhenoSS is an effective algorithm that makes disease prediction and performs patient clustering based on HPO concepts. PhenoSS uses the Gaussian copula technique by modeling the marginal prevalence of each HPO term for each disease and utilizes a multivariate normal distribution to link them together to account for term correlations. We utilized the OARD (open annotations for rare diseases) API for inferring the frequency of HPO terms in a diverse range of rare diseases. PhenoSS can calculate the phenotype similarity between any two patients for finding similar patients or clustering purposes, or between one patient and any candidate diseases for diagnosis support. 

The toolkit is implemented in Python. 

# Installation

### 1) Ensure to download several databases for this application:

**DIshIN** (required for ssmpy)
```
curl -L -O http://labs.rd.ciencias.ulisboa.pt/dishin/hp202506.db.gz
gunzip -N hp202506.db.gz
```
(you will find "hp.db" on the local directory)

**HPO-Disease Frequency**
- Access https://hpo.jax.org/data/annotations
- Download "GENES TO PHENOTYPE" and save in ./doc/database/

**OARD**
- There is no need to download anything here as the database will be called out during running. (https://rare.cohd.io/)

**MONDO**
- You may need this data if you need to convert from OMIM to MONDO or Orphanet to MONDO. Howewever, we directly provide the conversion files (./doc/database/omim_conversion.json | orphanet_conversion.json) for you (note that this can be outdated).
- Alternatively, download the mondo-edges.tsv for conversion.

### 2) Code environment:
Ensure that you install required packages:
```
pip install pandas numpy requests ssmpy scipy
```

### 3) You need to run 'python ./scripts/processing_hpo_frequency.py' to obtain ./doc/database/hpo_frequency.csv (one-time only)



# Tutorials

## Patient Clustering
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
The first input argument is the input file that contains the patient IDs and the HPO terms.
```
python ./scripts/similarity_score.py -input_dir [YOUR INPUT DIRECTORY] -output_dir [YOUR OUTPUT DIRECTORY]
```

The outputs of the argument can be found in the file 1_10_sim.
| pat1 | pat2 | hpo1 | hpo2 | similarity |
|------|------|------|------|-----------|
| 0_10 | 1_10 | HP_0004370;HP_0000280;HP_0002835;... | HP_0000483;HP_0002307;HP_0001090;... | 3.74 |
| 0_10 | 2_10 | HP_0004370;HP_0000280;HP_0002835;... | HP_0002194;HP_0001263;HP_0001684 | 8.09 |
| 1_10 | 2_10 | HP_0000483;HP_0002307;HP_0001090;... | HP_0002194;HP_0001263;HP_0001684;... | 2.45 |

## Disease prediction
PhenoSS extracts the diseases/phenotype frequencies from the Open Annotations for Rare Diseases (OARD) and Human Phenotype Ontology Databases. It takes in HPO terms of a list of patients and outputs the ranks of possible underlying diseases. 

Below is a sample input file:

```
P1	HP_0012759;HP_0000750;HP_0100022;HP_0000707;
P2	HP_0001270;HP_0012758;HP_0002066;HP_0011443;
P3	HP_0012758;HP_0002167;HP_0012638;HP_0000707;
```
You can use "PhenoSS_Codebook.ipynb" if you prefer the interactive browser. Otherwise to run PhenoSS, use the following command:

```
bash run_phenoss.sh \
  --inputfile data/patient_hpos.tsv \
  --outputfile results/phenoss_output.tsv \
  --mode hpo_first \
  --freq_assignment extrinsic_ic \
  --method Resnik \
  --hp_db_sqlite hp.db \
  --hpo_db_path ./databases/hpo_frequency.csv \
  --url https://rare.cohd.io/api \
  --dataset_id 2 \
  --gene_conversion \
  --gene_of_interest GATA2 \
  --gene_outfile results/gene_results.tsv
```
| Argument             | Required  | Default                    | Description                                                                                   |
| -------------------- | --------- | -------------------------- | --------------------------------------------------------------------------------------------- |
| `--inputfile`        | **Yes**   | —                          | Input file containing patient HPO phenotypes                                                  |
| `--outputfile`       | **Yes**   | —                          | Output ranking file                                                                           |
| `--mode`             | No        | `hpo_first`               | Candidate disease selection strategy (`oard_only`, `oard_first`, `hpodb_first`, `hpodb_only`) |
| `--freq_assignment`  | No        | `assumption`             | Frequency assignment method for HPO terms                                                     |
| `--method`           | No        | `Resnik`                   | Semantic similarity method                                                                    |
| `--hp_db_sqlite`     | No        | `hp.db`                    | SQLite database for HPO ontology                                                              |
| `--hpo_db_path`      | No        | `./databases/hpo_frequency.csv`        | HPO frequency table                                                                           |
| `--gene_conversion`  | No (flag) | Off                        | Convert diseases to genes in output. Note that diseases with unknown genes will be removed.                                                           |
| `--hpo_removal`  | No (flag) | Off                        | Remove less informative HPO terms such as HP:0000118 to increase the precision.                                                           |                                                     |
| `--limit`              | No        | `0` | Maximum number of diseases can be considered when predicting per sample (a list of HPO terms may be linked to thousands of diseases). Smaller is faster. Set default to be 0 means no limit.                                                                             |
| `--url`              | No        | `https://rare.cohd.io/api` | COHD API endpoint                                                                             |
| `--dataset_id`       | No        | `2`                        | Dataset ID for COHD                                                                           |
| `--gene_of_interest` | No        | empty                      | Specific gene to evaluate                                                                     |
| `--gene_outfile`     | No        | empty                      | Output file for gene results                                                                  |


The results consist of a list of MONDO diseases and the rankings and will be stored in 'outputFile' specified by the user. 

### Output Example
| patient_id  | disease_mondo | gene  | disease_id | score              | rank | 
| ----------- | ------------- | ----- | ---------- | ------------------ | ---- |
| PATIENT_001 | MONDO:0001234 | GENE1 | 80001234   | -0.941566671739919 | 1    | 
| PATIENT_001 | MONDO:0001234 | GENE2 | 80001234   | -0.941566671739919 | 1    | 
| PATIENT_001 | MONDO:0005678 | NA | 80005678   | -0.941503249828281 | 3    | 


## License

PhenoSS is distributed under the [MIT License by Wang Genomics Lab](https://wglab.mit-license.org/).
