import pandas as pd
import numpy as np
import json
with open('./databases/omim_conversion.json', 'r') as f:
    omim_dict = json.load(f)
with open('./databases/orphanet_conversion.json', 'r') as f:
    orphanet_dict = json.load(f)
def convert_frequency(freq):
    if pd.isna(freq):
        return 0.01 # very rare
    elif "HP" in freq:
        transform = {'HP:0040283':0.33,
        'HP:0040281':0.9,
        'HP:0040282':0.5,
        'HP:0040284':0.01,
        'HP:0040280':1}
        try:
            return transform[freq]
        except:
            return 0.01
    elif '/' in freq:
        slash_idx = freq.find('/')
        freq_num = float(freq[:slash_idx]) / float(freq[slash_idx+1:])
        return freq_num
    elif '%' in freq:
        slash_idx = freq.find('%')
        freq_num = float(freq[:slash_idx]) / 100
        return freq_num
    else:
        return 0.01 # very rare
def mondo_convert(x):
    if "OMIM" in x:
        if x in omim_dict:
            return ";".join(omim_dict[x])
        else:
            return np.nan
    elif "ORPHA" in x:
        if x in orphanet_dict:
            return ";".join(orphanet_dict[x])
        else:
            return np.nan
    else:
        np.nan

gene2phen = pd.read_csv('./databases/genes_to_phenotype.txt', sep = '\t')
gene2phen['frequency'] = gene2phen['frequency'].apply(convert_frequency)
gene2phen['mondo_id'] = gene2phen['disease_id'].apply(mondo_convert)
gene2phen["mondo_id"] = gene2phen["mondo_id"].str.split(";")
gene2phen = gene2phen.explode("mondo_id").reset_index(drop=True)
gene2phen = gene2phen[pd.notnull(gene2phen['mondo_id'])].reset_index(drop=True)
gene2phen.to_csv('./databases/hpo_frequency.csv', sep = '\t', index=False)