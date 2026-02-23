from __future__ import print_function # load print function in python3
from collections import defaultdict
import numpy as np
import pandas as pd
import ssmpy, sys, json, argparse, os
from tqdm import tqdm
ssmpy.ssm.mica = True
ssmpy.ssm.intrinsic = True
from copy import deepcopy
ssmpy.semantic_base("hp.db")
def shard_dict(merged_output, index, num_shards=30):
    """
    Split a dictionary into `num_shards` parts using stable modulo sharding.
    If index == num_shards, return the remaining (unassigned) keys.
    """
    assert 0 <= index < num_shards

    keys = sorted(merged_output.keys())  # deterministic order

    shard_keys = [k for i, k in enumerate(keys) if i % num_shards == index]

    return {k: merged_output[k] for k in shard_keys}
def auto_dict():
    return defaultdict(auto_dict)


# cache: HPO string -> ssmpy internal ID
hpo_id_cache = {}
def get_eid(hp):
    if hp not in hpo_id_cache:
        hpo_id_cache[hp] = ssmpy.get_id(hp)
    return hpo_id_cache[hp]

# cache: (eid1, eid2) -> Resnik score
resnik_cache = {}

def resnik(e1, e2):
    key = (e1, e2) if e1 <= e2 else (e2, e1)
    if key not in resnik_cache:
        resnik_cache[key] = ssmpy.ssm_resnik(e1, e2)
    return resnik_cache[key]
def bma_similarity(hpo_set1, hpo_set2):
    """
    Bidirectional Best Match Average using Resnik similarity
    """
    if len(hpo_set1) == 0 or len(hpo_set2) == 0:
        return np.nan

    # A -> B
    score1 = np.mean([
        max(resnik(get_eid(h1), get_eid(h2)) for h2 in hpo_set2)
        for h1 in hpo_set1
    ])

    # B -> A
    score2 = np.mean([
        max(resnik(get_eid(h2), get_eid(h1)) for h1 in hpo_set1)
        for h2 in hpo_set2
    ])

    return (score1 + score2) / 2

def main():
    parser = argparse.ArgumentParser(description="PhenoGPT2 Phenotype Recognizer and Normalizer")
    parser.add_argument("-index", "--index", type=int, help="Index identifier for saving")
    args = parser.parse_args()
    print("Read the data", flush=True)
    with open('/home/nguyenqm/projects/KnowledgeGraph/data/OMIM_diseases_KG.json', 'r') as f:
        all_kg = json.load(f)

    gmdb_data = pd.read_csv('/home/nguyenqm/projects/MedicalDBs/GMDB_v1.1.0/all_train_both.csv', sep = '\t')
    gmdb_data['OMIM'] = gmdb_data['OMIM'].apply(lambda x: str(int(x)) if pd.notnull(x) else x)
    gmdb_data = gmdb_data[pd.notnull(gmdb_data['OMIM'])]
    all_omims = gmdb_data['OMIM'].unique()
    df = pd.read_csv(
        "/home/nguyenqm/projects/MedicalDBs/OMIM/mimTitles.txt",
        skiprows=2,   # skip first 2 rows
        comment=None, # do NOT treat # as comment globally
        sep = '\t'
    )
    # Remove leading '#' and whitespace from column names
    df.columns = df.columns.str.lstrip("# ").str.strip()
    df = df[:-13]
    df['MIM Number'] = df['MIM Number'].apply(lambda x: str(int(x)))
    #remaining = list(set(all_omims) - set(all_kg.keys()))
    common_omims = list(set(all_omims) & set(all_kg.keys()))
    print(f"Number of OMIM Diseases: {len(common_omims)}")
    gmdb_dx = df[df['MIM Number'].isin(common_omims)].set_index('MIM Number')['Preferred Title; symbol'].to_dict()
    facial_organs = ['face', 'eye', 'neck', 'ear', 'mouth', 'eye', 
                    'nose','head', 'hair', 'dental', 'skull', 'cheek','forehead','lip']
    dis1 = auto_dict()
    #dis2 = auto_dict()
    print("Start getting data")
    comparedone = auto_dict()
    for omim,data in all_kg.items():
        #if omim in list(top10_omim.keys()):
        if omim in list(gmdb_dx.keys()):
            all_hpos = []
            for phen, phen_dict in data['phenotypes'].items():
                if phen_dict['polarity'] == 'present' and "HP" in phen_dict['hpo']:# and (any([x in phen_dict['organ'] for x in facial_organs])):
                    all_hpos.append(phen_dict['hpo'].replace(":","_"))
            dis1[omim] = list(set(all_hpos))
            #dis2[omim] = list(set(all_hpos))
    dis1 = dict(sorted(dis1.items()))
    dis2 = deepcopy(dis1)
    all_pats = sorted(dis1.keys())
    pat2idx = {p: i for i, p in enumerate(all_pats)}
    #remaining_omims = ['601471', '217980', '209850', '614066', '613970', '614202', '614104', '136500', '236100', '300082', '601853', '614098', '157900', '614080', '241550']
    # dis1 = {k:v for k,v in dis1.items() if k in remaining_omims}
    fileindex = args.index
    rows = []
    patients = sorted(dis1.keys())
    output_dir = '/home/nguyenqm/projects/github/PhenoSS/gmdb_similarity_all_hpos'
    os.makedirs(output_dir, exist_ok = True)
    out_file = f"{output_dir}/phenoss{fileindex}.csv"
    print("Start working", flush=True)
    for i, pat1 in enumerate(tqdm(patients, desc="pat1")):
        print(pat1, flush=True)
        i = pat2idx[pat1]
        hpo1 = dis1[pat1]#.split(";")

        for pat2 in tqdm(all_pats[i+1:]):  # ensures no (pat2, pat1)
        # for pat2 in tqdm(all_pats):  # ensures no (pat2, pat1)
            if pat1 == pat2:
                continue
            hpo2 = dis2[pat2]#.split(";")

            sim = bma_similarity(hpo1, hpo2)

            rows.append({
                "omim1": pat1,
                "omim2": pat2,
                "hpo1": ";".join(hpo1),
                "hpo2": ";".join(hpo2),
                "similarity": sim
            })

    df = pd.DataFrame(rows)
    df.to_csv(out_file, sep = '\t', index=False)
    
if __name__ == "__main__":
    main()