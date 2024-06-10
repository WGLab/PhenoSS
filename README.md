# PhenoSS: Phenotype semantic similarity-based approach for rare disease prediction and patient clustering 

## Introduction
PhenoSS is an effective algorithm that makes disease prediction and performs patient clustering based on HPO concepts. PhenoSS uses the Gaussian copula technique by modeling the marginal prevalence of each HPO term for each disease and utilizes a multivariate normal distribution to link them together to account for term correlations. We utilized the OARD (open annotations for rare diseases) API for inferring the frequency of HPO terms in a diverse range of rare diseases. PhenoSS can calculate the phenotype similarity between any two patients for finding similar patients or clustering purposes, or between one patient and any candidate diseases for diagnosis support. 

The tookit is implemented in Python. 

## Installation

## Datasets
#### Human Phenotype Ontology (HPO): 
https://hpo.jax.org/
#### Open Annotations for Rare Diseases (OARD):
https://rare.cohd.io/
