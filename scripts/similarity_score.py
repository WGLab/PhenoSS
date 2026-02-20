from __future__ import print_function  # load print function in python3
from collections import defaultdict
import numpy as np
import pandas as pd
import ssmpy, sys, json, argparse, os
from tqdm import tqdm
from copy import deepcopy
# ----------------------------
# PhenoSS configuration
# ----------------------------
ssmpy.ssm.mica = True
ssmpy.ssm.intrinsic = True
ssmpy.semantic_base("hp.db")

# ----------------------------
# Utility functions (UNCHANGED)
# ----------------------------
def shard_dict(merged_output, index, num_shards=30):
    """
    Split a dictionary into `num_shards` parts using stable modulo sharding.
    """
    assert 0 <= index < num_shards
    keys = sorted(merged_output.keys())
    shard_keys = [k for i, k in enumerate(keys) if i % num_shards == index]
    return {k: merged_output[k] for k in shard_keys}

def auto_dict():
    return defaultdict(auto_dict)

# ----------------------------
# Caching (UNCHANGED)
# ----------------------------
hpo_id_cache = {}
def get_eid(hp):
    if hp not in hpo_id_cache:
        hpo_id_cache[hp] = ssmpy.get_id(hp)
    return hpo_id_cache[hp]

resnik_cache = {}
def resnik(e1, e2):
    key = (e1, e2) if e1 <= e2 else (e2, e1)
    if key not in resnik_cache:
        resnik_cache[key] = ssmpy.ssm_resnik(e1, e2)
    return resnik_cache[key]

# ==========================================================
# ✅ PhenoSS similarity (matches original implementation)
# ==========================================================
def phenoss_directional(hpo_src, hpo_tgt):
    """
    Exact PhenoSS directional similarity:
    diff = (sum_{t in src} max_{t' in tgt} Resnik(t, t')) / |src|

    This matches the original code:
        totaldiff += maxdiff
        diff = totaldiff / numcompare
    """
    if len(hpo_src) == 0 or len(hpo_tgt) == 0:
        return np.nan

    totaldiff = 0.0
    numcompare = 0

    for h1 in hpo_src:
        e1 = get_eid(h1)

        # same as building alldiff then np.amax(alldiff)
        maxdiff = max(resnik(e1, get_eid(h2)) for h2 in hpo_tgt)

        totaldiff += maxdiff
        numcompare += 1

    return totaldiff / numcompare if numcompare > 0 else np.nan


def phenoss_similarity(hpo_set1, hpo_set2, identity=False):
    """
    Symmetric PhenoSS (as in author's Perl script):
        sim = 0.5 * (sim_{P1->P2} + sim_{P2->P1})
    """
    if len(hpo_set1) == 0 or len(hpo_set2) == 0:
        return np.nan
    if identity:
        sim_p1_to_p2 = phenoss_directional(hpo_set1, hpo_set2)
        return sim_p1_to_p2
    sim_p1_to_p2 = phenoss_directional(hpo_set1, hpo_set2)
    sim_p2_to_p1 = phenoss_directional(hpo_set2, hpo_set1)

    return 0.5 * (sim_p1_to_p2 + sim_p2_to_p1)

# ----------------------------
# Main (DATA PIPELINE UNCHANGED)
# ----------------------------
def main():
    parser = argparse.ArgumentParser(description="PhenoSS patient similarity")
    parser.add_argument("--input_file", required=True, help="Input TSV file (patient_id  HPOs)")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Reading input file...", flush=True)

    # --------------------------------------------------
    # READ INPUT
    # --------------------------------------------------
    df = pd.read_csv(args.input_file, sep="\t", header=None, names=["patient", "hpos"])

    # convert to dictionary
    pat2hpo = {}
    for _, row in df.iterrows():
        hpos = [h.strip() for h in str(row["hpos"]).split(";") if h.strip() != ""]
        pat2hpo[row["patient"]] = hpos

    patients = list(pat2hpo.keys())

    print(f"Loaded {len(patients)} patients", flush=True)

    # --------------------------------------------------
    # COMPUTE PAIRWISE SIMILARITY
    # --------------------------------------------------
    rows = []

    for i, pat1 in enumerate(tqdm(patients, desc="pat1")):
        hpo1 = pat2hpo[pat1]

        for pat2 in patients[i+1:]:
            hpo2 = pat2hpo[pat2]

            sim = phenoss_similarity(hpo1, hpo2)

            rows.append({
                "pat1": pat1,
                "pat2": pat2,
                "hpo1": ";".join(hpo1),
                "hpo2": ";".join(hpo2),
                "similarity": sim
            })

    # --------------------------------------------------
    # SAVE OUTPUT
    # --------------------------------------------------
    df_out = pd.DataFrame(rows)
    out_file = os.path.join(args.output_dir, "phenoss_similarity.tsv")
    df_out.to_csv(out_file, sep="\t", index=False)

    print(f"Saved: {out_file}", flush=True)


if __name__ == "__main__":
    main()
