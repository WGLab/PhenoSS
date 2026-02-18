#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reorganized + argparse wrapper for the provided script.

Modes:
  a) oard_only      : only use OARD (no HPO_db updates; no HPO_db diseases; no HPO_db pair frequencies)
  b) oard_first     : prioritize OARD over HPO_db (HPO update + pair freq preference)
  c) hpodb_first    : prioritize HPO_db over OARD
  d) hpodb_only     : only use HPO_db (no OARD calls)

freq_assignment (used when HPO_db is enabled in mode b/c/d):
  mean, median, max,
  extrinsic_ic, intrinsic_ic, ic (ssmpy.ssm.information_content),
  assumption (fallback/else branch: 0.01)
"""

import argparse
from collections import defaultdict, Counter
from tqdm import tqdm
import numpy as np
import pandas as pd
import requests
import ssmpy
from scipy.stats import multivariate_normal, norm
import warnings
warnings.filterwarnings("ignore")

MVN_CACHE = {}
# =========================================================
# SSMPY SETUP
# =========================================================
def init_ssmpy(hp_db_sqlite: str) -> None:
    ssmpy.ssm.mica = True
    ssmpy.ssm.intrinsic = True
    ssmpy.semantic_base(hp_db_sqlite)


# =========================================================
# LOAD DATA
# =========================================================
def load_oard(url: str, dataset_id: str = "2"):
    # HPO domain
    params = {"dataset_id": dataset_id, "domain_id": "phenotypes"}
    oard_data = pd.DataFrame(
        requests.get(url + "/frequencies/mostFrequency", params=params, verify=False)
        .json()["results"]
    )

    hpo_oard = list(set(oard_data["concept_code"]))
    oard_hpo_freq_map = oard_data.set_index("concept_code")["concept_frequency"].to_dict()
    oard_hpo_code2id = oard_data.set_index("concept_code")["concept_id"].to_dict()

    # Disease domain
    params = {"dataset_id": dataset_id, "domain_id": "diseases"}
    domain_disease_df = pd.DataFrame(
        requests.get(url + "/frequencies/mostFrequency", params=params, verify=False)
        .json()["results"]
    )

    oard_disease = list(set(domain_disease_df["concept_id"]))

    return {
        "oard_data": oard_data,
        "hpo_oard": hpo_oard,
        "oard_hpo_freq_map": oard_hpo_freq_map,
        "oard_hpo_code2id": oard_hpo_code2id,
        "domain_disease_df": domain_disease_df,
        "oard_disease": oard_disease
    }


def load_hpo_db(hpo_db_path: str):
    hpo_db = pd.read_csv(hpo_db_path, sep="\t")
    hpo_db = hpo_db.rename(
        columns={
            "hpo_id": "concept_code_1",
            "mondo_id": "concept_code_2",
            "frequency": "concept_frequency",
        }
    )
    hpo_db["concept_id_1"] = hpo_db["concept_code_1"].apply(lambda x: int(x.replace("HP:", "9")))
    hpo_db["concept_id_2"] = hpo_db["concept_code_2"].apply(lambda x: int(x.replace("MONDO:", "8")))

    hpo_db_freq_lookup = dict(
        zip(
            zip(hpo_db["concept_code_1"], hpo_db["concept_id_2"]),
            hpo_db["concept_frequency"],
        )
    )
    hpo_db_index = hpo_db.groupby("concept_code_1")["concept_id_2"].apply(list).to_dict()

    return {
        "hpo_db": hpo_db,
        "hpo_db_freq_lookup": hpo_db_freq_lookup,
        "hpo_db_index": hpo_db_index,
    }


# =========================================================
# MARGINAL FREQ FROM HPO_DB
# =========================================================
def build_hpo_db_marginal_freq(hpo_db: pd.DataFrame, freq_assignment: str):
    """
    Returns a function: get_freq(hpo_code) -> float
    Computes marginal frequency only for requested HPOs (lazy evaluation).
    """

    hpo_db_group = hpo_db.groupby("concept_code_1")["concept_frequency"]
    cache = {}

    def get_freq(hpo_code: str):

        if hpo_code in cache:
            return cache[hpo_code]

        if freq_assignment == "mean":
            val = float(hpo_db_group.mean().get(hpo_code, 0.01))

        elif freq_assignment == "median":
            val = float(hpo_db_group.median().get(hpo_code, 0.01))

        elif freq_assignment == "max":
            val = float(hpo_db_group.max().get(hpo_code, 0.01))

        elif freq_assignment in {"extrinsic_ic", "intrinsic_ic", "ic"}:
            hid = ssmpy.get_id(hpo_code.replace(":", "_"))
            if hid != -1:
                if freq_assignment == "extrinsic_ic":
                    ic = ssmpy.ssm.information_content_extrinsic(hid)
                elif freq_assignment == "intrinsic_ic":
                    ic = ssmpy.ssm.information_content_intrinsic(hid)
                else:
                    ic = ssmpy.ssm.information_content(hid)
                val = float(np.exp(-ic))
            else:
                val = 0.01

        else:  # assumption
            val = 0.01

        cache[hpo_code] = val
        return val

    return get_freq
def get_hpo_ic_weight(hpo_code, caches, ic_type="extrinsic_ic"):

    if not hasattr(caches, "hpo_ic_cache"):
        caches.hpo_ic_cache = {}

    key = (hpo_code, ic_type)
    if key in caches.hpo_ic_cache:
        return caches.hpo_ic_cache[key]

    hid = ssmpy.get_id(hpo_code.replace(":", "_"))
    if hid == -1:
        w = 0.0
    else:
        if ic_type == "extrinsic_ic":
            w = float(ssmpy.ssm.information_content_extrinsic(hid))
        elif ic_type == "intrinsic_ic":
            w = float(ssmpy.ssm.information_content_intrinsic(hid))
        else:
            w = float(ssmpy.ssm.information_content(hid))

    caches.hpo_ic_cache[key] = w
    return w


# =========================================================
# CACHES + SIMILARITY
# =========================================================
class Caches:
    def __init__(self):
        self.oard_hpo_api_cache = {}
        self.hpo_freq_cache = {}
        self.ssmpy_id_cache = {}

    def get_ssmpy_id(self, h: str) -> int:
        if h not in self.ssmpy_id_cache:
            self.ssmpy_id_cache[h] = ssmpy.get_id(h)
        return self.ssmpy_id_cache[h]


def calc_sim(hpo1: str, hpo2: str, caches: Caches, method: str = "Resnik"):
    e1 = caches.get_ssmpy_id(hpo1)
    e2 = caches.get_ssmpy_id(hpo2)
    if e1 == -1 or e2 == -1:
        return 0
    elif method == "Resnik":
        return ssmpy.ssm_resnik(e1, e2)
    elif method == "Lin":
        return ssmpy.ssm_lin(e1, e2)
    elif method == "JC":
        return ssmpy.ssm_jiang_conrath(e1, e2)
    elif method == "IC":
        return ssmpy.ssm_lin(e1, e2) * (1 - 1 / (1 + ssmpy.ssm_resnik(e1, e2)))


# =========================================================
# ODDS + UPDATE_HPO (mode-dependent priority)
# =========================================================
def calc_odd_oard_fast(disease_idx, hpos, hpos_freq, rho, map2freq, x_base, cov_mat, mvn=None):
    hpos = [each.replace("_", ":") for each in hpos]
    p = 10 ^ -6
    odd_disease = p / (1 - p)

    x_disease = x_base.copy()
    x_noDisease = x_base.copy()

    for i, hpo in enumerate(hpos):
        p_yes = map2freq.get((hpo.replace(":", "_"), disease_idx), -1)
        if p_yes != -1:
            x_disease[i] = norm.ppf(p_yes)
            p_no = (hpos_freq[i] - p_yes * p) / (1 - p)
            if p_no < 0:
                p_no = 0
            x_noDisease[i] = norm.ppf(p_no)

    if mvn is None:
        joint_prob1 = multivariate_normal.cdf(x_disease, cov=cov_mat)
        joint_prob2 = multivariate_normal.cdf(x_noDisease, cov=cov_mat)
    else:
        joint_prob1 = mvn.cdf(x_disease)
        joint_prob2 = mvn.cdf(x_noDisease)


    return odd_disease * joint_prob1 / joint_prob2


def update_hpo(
    hpos,
    *,
    method: str,
    mode: str,
    caches: Caches,
    # OARD
    hpo_oard=None,
    oard_hpo_freq_map=None,
    # HPO_db marginal
    hpo_db_get_marginal_freq=None,
    hpo_default_frequency: float = 0.01,
):
    """
    Returns: (new_hpos, hpo_freq)
    Priority is controlled by mode:
      - oard_only:     OARD -> map-to-nearest-OARD
      - hpodb_only:    HPO_db -> else default
      - oard_first:    OARD -> HPO_db -> map-to-nearest-OARD
      - hpodb_first:   HPO_db -> OARD -> map-to-nearest-OARD
    """
    new_hpos = []
    hpo_freq = []

    for hpo in tqdm(hpos, position=1,  dynamic_ncols=True, leave=True):
        hpo_code = hpo.replace("_", ":")

        if mode == "oard_only":
            if hpo_code in (oard_hpo_freq_map or {}):
                new_hpos.append(hpo)
                hpo_freq.append(float(oard_hpo_freq_map[hpo_code]))
                continue
            # map to nearest OARD (Resnik)
            if hpo not in caches.hpo_freq_cache:
                sim_scores = [calc_sim(each.replace(":", "_"), hpo, caches, method) for each in (hpo_oard or [])]
                index_max = int(np.argmax(sim_scores)) if len(sim_scores) else 0
                mapped_hpo = hpo_oard[index_max].replace(":", "_")
                caches.hpo_freq_cache[hpo] = mapped_hpo
            else:
                mapped_hpo = caches.hpo_freq_cache[hpo]

            new_hpos.append(mapped_hpo)
            mapped_code = mapped_hpo.replace("_", ":")
            hpo_freq.append(float((oard_hpo_freq_map or {}).get(mapped_code, hpo_default_frequency)))
            continue

        if mode == "hpodb_only":
            if hpo_db_get_marginal_freq is not None:
                val = hpo_db_get_marginal_freq(hpo_code)
                if val is not None:
                    new_hpos.append(hpo)
                    hpo_freq.append(float(val))
                    continue
            else:
                new_hpos.append(hpo)
                hpo_freq.append(float(hpo_default_frequency))
            continue

        if mode == "oard_first":
            # 1) OARD
            if hpo_code in (oard_hpo_freq_map or {}):
                new_hpos.append(hpo)
                hpo_freq.append(float(oard_hpo_freq_map[hpo_code]))
                continue
            # 2) HPO_db
            if hpo_db_get_marginal_freq is not None:
                val = hpo_db_get_marginal_freq(hpo_code)
                if val is not None:
                    new_hpos.append(hpo)
                    hpo_freq.append(float(val))
                    continue

            # 3) map to nearest OARD
            if hpo not in caches.hpo_freq_cache:
                sim_scores = [calc_sim(each.replace(":", "_"), hpo, caches, method) for each in (hpo_oard or [])]
                index_max = int(np.argmax(sim_scores)) if len(sim_scores) else 0
                mapped_hpo = hpo_oard[index_max].replace(":", "_")
                caches.hpo_freq_cache[hpo] = mapped_hpo
            else:
                mapped_hpo = caches.hpo_freq_cache[hpo]

            new_hpos.append(mapped_hpo)
            mapped_code = mapped_hpo.replace("_", ":")
            hpo_freq.append(float((oard_hpo_freq_map or {}).get(mapped_code, hpo_default_frequency)))
            continue

        if mode == "hpodb_first":
            # 1) HPO_db
            if hpo_db_get_marginal_freq is not None:
                val = hpo_db_get_marginal_freq(hpo_code)
                if val is not None:
                    new_hpos.append(hpo)
                    hpo_freq.append(float(val))
                    continue
            # 2) OARD
            if hpo_code in (oard_hpo_freq_map or {}):
                new_hpos.append(hpo)
                hpo_freq.append(float(oard_hpo_freq_map[hpo_code]))
                continue
            # 3) map to nearest OARD
            if hpo not in caches.hpo_freq_cache:
                sim_scores = [calc_sim(each.replace(":", "_"), hpo, caches, method) for each in (hpo_oard or [])]
                index_max = int(np.argmax(sim_scores)) if len(sim_scores) else 0
                mapped_hpo = hpo_oard[index_max].replace(":", "_")
                caches.hpo_freq_cache[hpo] = mapped_hpo
            else:
                mapped_hpo = caches.hpo_freq_cache[hpo]

            new_hpos.append(mapped_hpo)
            mapped_code = mapped_hpo.replace("_", ":")
            hpo_freq.append(float((oard_hpo_freq_map or {}).get(mapped_code, hpo_default_frequency)))
            continue

        # Fallback (should not happen)
        new_hpos.append(hpo)
        hpo_freq.append(float(hpo_default_frequency))

    return new_hpos, hpo_freq


# =========================================================
# PATIENT INPUT
# =========================================================
def load_patients(inputfile: str, *, update_fn):
    hpo_patients = {}
    freq_patients = {}
    with open(inputfile) as infile:
        for line in tqdm(infile, position=0, desc = 'Entering patient info'):
            parts = line.strip().split()
            pid = parts[0]
            hpos = parts[1].split(";")[:-1]
            new_hpos, new_freq = update_fn(hpos)
            hpo_patients[pid] = new_hpos
            freq_patients[pid] = new_freq
    return hpo_patients, freq_patients


# =========================================================
# BUILD disease2freq + hpo2diseases (mode-dependent)
# =========================================================
def build_disease_maps(
    *,
    mode: str,
    caches: Caches,
    url: str,
    dataset_id: str,
    # patient data
    hpo_patients: dict,
    # OARD
    oard_hpo_code2id=None,
    # HPO_db
    hpo_db_index=None,
    hpo_db_freq_lookup=None,
    # preference
    prefer: str,  # "oard" or "hpodb"
):
    disease2freq = {}
    hpo2diseases = defaultdict(set)

    final_diseases = set()
    all_patient_hpos = set([h for lst in hpo_patients.values() for h in lst])

    for each_hpo in tqdm(all_patient_hpos, desc='Retrieving Frequency'):
        hpo_key = each_hpo.replace(":", "_")
        each_hpo_code = each_hpo.replace("_", ":")

        # ---- collect candidate diseases ----
        if mode in {"hpodb_only", "hpodb_first", "oard_first"}:
            final_diseases.update((hpo_db_index or {}).get(each_hpo_code, []))

        # ---- OARD pull for this HPO ----
        oard_lookup = {}
        oard_final_dict = set()
        hpo_id = None

        if mode in {"oard_only", "oard_first", "hpodb_first"}:
            hpo_id = (oard_hpo_code2id or {}).get(each_hpo_code, None)
            if hpo_id:
                if hpo_id not in caches.oard_hpo_api_cache:
                    params = {"dataset_id": dataset_id, "concept_id": hpo_id}
                    result_df = pd.DataFrame(
                        requests.get(url + "/frequencies/mostFrequency", params=params, verify=False)
                        .json()["results"]
                    )
                    caches.oard_hpo_api_cache[hpo_id] = result_df
                else:
                    result_df = caches.oard_hpo_api_cache[hpo_id]

                if len(result_df) > 0:
                    result_df = result_df[result_df["concept_id_2"].astype(str).str.startswith("8")]
                    oard_lookup = dict(
                        zip(
                            zip(result_df["concept_id_1"], result_df["concept_id_2"]),
                            result_df["concept_frequency"],
                        )
                    )
                    oard_final_dict = set([d for (_, d) in oard_lookup.keys()])

        if mode in {"oard_only", "oard_first", "hpodb_first"}:
            final_diseases.update(oard_final_dict)

        # ---- fill disease2freq + hpo2diseases ----
        for dis in final_diseases:
            # preference for pair frequency
            if mode == "oard_only":
                val = oard_lookup.get((hpo_id, dis), -1) if hpo_id else -1

            elif mode == "hpodb_only":
                val = (hpo_db_freq_lookup or {}).get((each_hpo_code, dis), -1)

            else:
                # mixed
                if prefer == "oard":
                    val = oard_lookup.get((hpo_id, dis), -1) if hpo_id else -1
                    if val == -1:
                        val = (hpo_db_freq_lookup or {}).get((each_hpo_code, dis), -1)
                else:
                    val = (hpo_db_freq_lookup or {}).get((each_hpo_code, dis), -1)
                    if val == -1:
                        val = oard_lookup.get((hpo_id, dis), -1) if hpo_id else -1

            disease2freq[(hpo_key, dis)] = val
            hpo2diseases[hpo_key].add(dis)

    return disease2freq, hpo2diseases


# =========================================================
# RANK PATIENTS
# =========================================================
def get_mvn_for_dimension(n, cov_mat):
    if n not in MVN_CACHE:
        MVN_CACHE[n] = multivariate_normal(mean=np.zeros(n), cov=cov_mat)
    return MVN_CACHE[n]
def computing_phenoss(hpo_patients, freq_patients, disease2freq, hpo2diseases, args, caches):
    rank_patients = {}
    for pid, patient_hpos in tqdm(hpo_patients.items(), position=0, desc = 'Computing the score'):
        patient_freq = freq_patients[pid]

        x_base = np.array([norm.ppf(p) for p in patient_freq])
        n = len(patient_hpos)

        cov_mat = 0.01 * np.ones([n, n])
        np.fill_diagonal(cov_mat, 1)

        mvn = get_mvn_for_dimension(n, cov_mat)
        # ---------- IC-weighted candidate filtering ----------
        if args.mode == 'oard_only':
            candidate_diseases = set()
            for h in patient_hpos:
                candidate_diseases.update(hpo2diseases.get(h, []))
        else:
            if len(patient_hpos) <= 1:
                candidate_diseases = set()
                for h in patient_hpos:
                    candidate_diseases.update(hpo2diseases.get(h, []))
            else:

                disease_score = defaultdict(float)
                disease_count = defaultdict(int)

                for h in patient_hpos:
                    h_code = h.replace("_", ":")
                    w = get_hpo_ic_weight(h_code, caches)

                    for d in hpo2diseases.get(h, []):
                        disease_score[d] += w
                        disease_count[d] += 1

                # sort by IC-weighted score
                sorted_items = sorted(
                    disease_score.items(),
                    key=lambda x: (x[1], disease_count[x[0]]),
                    reverse=True
                )

                TARGET = 1500
                HARD_CAP = 3000

                if len(sorted_items) <= TARGET:
                    candidate_diseases = [d for d, _ in sorted_items]
                else:
                    cutoff_score = sorted_items[TARGET-1][1]

                    candidate_diseases = [
                        d for d, s in sorted_items
                        if s >= cutoff_score
                    ]

                    if len(candidate_diseases) > HARD_CAP:
                        candidate_diseases = [d for d, _ in sorted_items[:HARD_CAP]]
        scores = {}
        for dis in candidate_diseases:
            scores[dis] = calc_odd_oard_fast(
                dis,
                patient_hpos,
                patient_freq,
                0.01,
                disease2freq,
                x_base,
                cov_mat,
                mvn=mvn
            )

        rank_patients[pid] = scores
        print("finished patient:", pid)

    return rank_patients


# =========================================================
# OUTPUT
# =========================================================
def build_output_df(rank_patients, hpo_db: pd.DataFrame, gene_conversion: bool):
    mondo2gene = (
        hpo_db.groupby("concept_code_2")["gene_symbol"]
        .apply(lambda x: list(set(x.dropna())))
        .to_dict()
    )

    records = []
    for pid, dis2score in tqdm(rank_patients.items(), desc = 'Finalizing output'):
        for dis, score in dis2score.items():
            mondo_code = "MONDO:" + str(dis)[1:]
            records.append(
                {
                    "patient_id": pid,
                    "disease_mondo": mondo_code,
                    "gene": mondo2gene.get(mondo_code, np.nan),  # list
                    "disease_id": dis,
                    "score": score,
                }
            )

    df_out = pd.DataFrame(records)
    if gene_conversion:
        df_out = df_out.explode("gene", ignore_index=True)

        # remove NaN genes
        df_out = df_out[df_out["gene"].notna()].copy()

    # sort for consistency
    df_out.sort_values(["patient_id", "score"], ascending=[True, True], inplace=True, ignore_index=True)
    df_out["score"] = pd.to_numeric(df_out["score"], errors="coerce")

    # rerank per patient (same score -> same rank)
    df_out["rank"] = (
        df_out.groupby("patient_id")["score"].rank(method="min", ascending=True).astype(int)
    ).reset_index(drop=True)

    df_out["sample_size"] = len(df_out)
    return df_out


# =========================================================
# MAIN
# =========================================================
def parse_args():
    p = argparse.ArgumentParser(description="PhenoSS ranking runner (reorganized).")

    p.add_argument("--inputfile", required=True, type=str, help="Path to patient input file.")
    p.add_argument("--outputfile", required=True, type=str, help="Path to output TSV/CSV.")

    p.add_argument(
        "--mode",
        required=True,
        choices=["oard_only", "oard_first", "hpodb_first", "hpodb_only"],
        help=(
            "Evaluation mode: "
            "oard_only (a), oard_first (b), hpodb_first (c), hpodb_only (d)."
        ),
    )

    p.add_argument(
        "--freq_assignment",
        default="extrinsic_ic",
        choices=["mean", "median", "max", "extrinsic_ic", "intrinsic_ic", "ic", "assumption"],
        help="How to compute marginal HPO frequencies from HPO_db when HPO_db is enabled.",
    )

    p.add_argument("--method", default="Resnik", type=str, help="Similarity method for mapping to OARD terms.")
    p.add_argument("--hp_db_sqlite", default="hp.db", type=str, help="Path to hp.db for ssmpy.")
    p.add_argument("--hpo_db_path", default="hpo_frequency.csv", type=str, help="Path to hpo_frequency.csv (tab-separated).")
    p.add_argument("--gene_conversion", action="store_true", help="Convert Mondo IDs to Genes")

    p.add_argument("--url", default="https://rare.cohd.io/api", type=str, help="OARD API base url.")
    p.add_argument("--dataset_id", default="2", type=str, help="OARD dataset_id (string).")

    # Optional debug extraction like your last lines
    p.add_argument("--gene_of_interest", default=None, type=str, help="If set, save rows for this gene.")
    p.add_argument("--gene_outfile", default=None, type=str, help="If set, write the gene subset to this path.")

    return p.parse_args()

def main():
    args = parse_args()
    print("Initializing HPO Database")
    init_ssmpy(args.hp_db_sqlite)

    caches = Caches()

    # Load HPO_db (always loaded because output step needs mondo2gene from it in your code)
    print("Loading HPO Database")
    hpo_db_pack = load_hpo_db(args.hpo_db_path)
    hpo_db = hpo_db_pack["hpo_db"]
    hpo_db_freq_lookup = hpo_db_pack["hpo_db_freq_lookup"]
    hpo_db_index = hpo_db_pack["hpo_db_index"]

    # Load OARD only if mode needs it
    print("Loading OARD Database")
    if args.mode in {"oard_only", "oard_first", "hpodb_first"}:
        oard_pack = load_oard(args.url, dataset_id=args.dataset_id)
        hpo_oard = oard_pack["hpo_oard"]
        oard_hpo_freq_map = oard_pack["oard_hpo_freq_map"]
        oard_hpo_code2id = oard_pack["oard_hpo_code2id"]
    else:
        hpo_oard = []
        oard_hpo_freq_map = {}
        oard_hpo_code2id = {}

    # Build HPO_db marginal freq only if mode uses HPO_db during update_hpo
    if args.mode in {"hpodb_only", "hpodb_first", "oard_first"}:
        hpo_db_get_marginal_freq = build_hpo_db_marginal_freq(hpo_db, args.freq_assignment)
    else:
        hpo_db_get_marginal_freq = {}

    hpo_default_frequency = 0.01

    # update_hpo wrapper with mode + resources
    def _update_fn(hpos):
        return update_hpo(
            hpos,
            method=args.method,
            mode=args.mode,
            caches=caches,
            hpo_oard=hpo_oard,
            oard_hpo_freq_map=oard_hpo_freq_map,
            hpo_db_get_marginal_freq=hpo_db_get_marginal_freq,
            hpo_default_frequency=hpo_default_frequency,
        )

    # Load patients
    print("Loading Patient Information")
    hpo_patients, freq_patients = load_patients(args.inputfile, update_fn=_update_fn)

    # disease2freq / hpo2diseases
    if args.mode == "oard_first":
        prefer = "oard"
    elif args.mode == "hpodb_first" or args.mode == 'hpodb_only':
        prefer = "hpodb"
    else:
        prefer = "oard"

    print("Building HPO-Disease Mapping")
    disease2freq, hpo2diseases = build_disease_maps(
        mode=args.mode,
        caches=caches,
        url=args.url,
        dataset_id=args.dataset_id,
        hpo_patients=hpo_patients,
        oard_hpo_code2id=oard_hpo_code2id,
        hpo_db_index=hpo_db_index,
        hpo_db_freq_lookup=hpo_db_freq_lookup,
        prefer=prefer,
    )

    # Compute PhenoSS scores
    print("Computing PhenoSS scores")
    phenoss_computation = computing_phenoss(hpo_patients, freq_patients, disease2freq, hpo2diseases, args, caches)

    # Output dataframe
    print("Building output")
    df_out = build_output_df(phenoss_computation, hpo_db, args.gene_conversion)

    # Save
    df_out.to_csv(args.outputfile, sep="\t", index=False)
    print("Saved results to", args.outputfile)

    # Optional gene subset (keeps your “LMNA” style check without hard-coded paths)
    if args.gene_of_interest is not None:
        print("Found the gene of interest. Finding its ranking...")
        x = df_out[df_out["gene"] == args.gene_of_interest].reset_index()
        print(x.head())
        if args.gene_outfile is not None:
            x.to_csv(args.gene_outfile, sep="\t", index=False)
            print("Saved gene subset to", args.gene_outfile)


if __name__ == "__main__":
    main()
