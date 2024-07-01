


from scipy.stats import norm
from scipy.stats import multivariate_normal
import ssmpy, sys

ssmpy.ssm.mica = True
ssmpy.ssm.intrinsic = True

ssmpy.semantic_base("hp.db")

import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io
url = 'https://rare.cohd.io/api'

# get the list of all pheonotypes with frequency annotated in OARD
params = {"dataset_id" : "2", "domain_id": "phenotypes"}
response = requests.get(url + '/frequencies/mostFrequency',params= params, verify=False)
response_json = response.json()
domain_counts_df = pd.DataFrame(response_json['results'])
oard_data = domain_counts_df
hpo_oard = list(set(oard_data['concept_code']))

params = {"dataset_id" : "2", "domain_id": "diseases"}
response = requests.get(url + '/frequencies/mostFrequency',params= params, verify=False)
response_json = response.json()
domain_disease_df = pd.DataFrame(response_json['results'])
oard_disease = list(set(domain_disease_df['concept_id']))

freq_map = {'very rare': 0.01, 'rare': 0.05, 'occasional': 0.075, 'main': 0.25, 'frequent': 0.33, 'typical':0.5, 'common':0.5, 'hallmark':0.9}

def calc_sim(hpo1, hpo2, method):
    # calculate the similarity score Resnik, Lin, JC, relevence, IC, graph_IC
    e1 = ssmpy.get_id(hpo1)
    #print(e1)
    e2 = ssmpy.get_id(hpo2)
    #print(e2)
    if e1==-1 or e2==-1:
        return 0
    elif method == 'Resnik':
        return ssmpy.ssm_resnik(e1,e2)
    elif method == 'Lin':
        return ssmpy.ssm_lin(e1, e2)
    elif method == 'JC':
        return ssmpy.ssm_jiang_conrath(e1, e2)
    elif method == 'IC':
        return ssmpy.ssm_lin(e1, e2)*(1-1/(1+ssmpy.ssm_resnik(e1,e2)))
        


def give_freq(hpo, table):  
    # assume that hpo is in the table

    freq = table[table['Term ID']==hpo.replace("_", ":")]['Frequency'].to_string(index=False)
    
    if '\n' in freq: # in case there are two frequencies corresponding to the same hpo
        freq = freq[:freq.find("\n")]
    freq_num = 0

    if '/' in freq:
        slash_idx = freq.find('/')
        freq_num = float(freq[:slash_idx]) / float(freq[slash_idx+1:])
    elif '%' in freq:
        slash_idx = freq.find('%')
        freq_num = float(freq[:slash_idx]) / 100
    else:
        freq_num = freq_map[freq.lower()]
                    
    if freq_num == 0:
        print("Error: Failed to find frequency")

    return freq_num




# so here the input list of hpos assume that their frequencies are all annotated in OARD
def calc_odd_oard(disease_idx, hpos, hpos_freq, rho,map2freq):  # disease_idx from 1 to 44
    hpos = [each.replace("_", ":") for each in hpos]
    #p = disease_freq[disease_idx]
    p = 10^-6
    odd_disease = p/(1-p)
    x_disease = []
    x_noDisease = []

    # assume if hpo not contained in table, then P(hpo|disease) = P(hpo|no disease) = P(hpo)

    for i in range(len(hpos)):
        hpo = hpos[i]
        p_yes = map2freq[(hpo.replace(":", "_"),disease_idx)]
        if p_yes == -1:
            x_disease.append(norm.ppf(hpos_freq[i]))
            x_noDisease.append(norm.ppf(hpos_freq[i]))
            #print(hpo)
        else:
            x_disease.append(norm.ppf(p_yes))
            #print(hpo)
            # use total probability law to calculate P(hpo|no disease)
            p_no = (hpos_freq[i] - p_yes * p)/(1-p)
            if p_no <0: p_no = 0
            #print(i)
            #print(p_no)
            #print(norm.ppf(p_no))
            x_noDisease.append(norm.ppf(p_no))
            
    #print(x_disease)

    # form the exchangable structure for the correlation matrix
    cov_mat = rho * np.ones([len(hpos), len(hpos)])
    n = cov_mat.shape[0]
    cov_mat[range(n), range(n)] = 1

    # multivariable gaussian 
    #print(x_disease)
    #print(x_noDisease)
    joint_prob1 = multivariate_normal.cdf(x_disease, cov = cov_mat)
    joint_prob2 = multivariate_normal.cdf(x_noDisease, cov = cov_mat)

    #print(joint_prob1)
    #print(joint_prob2)

    return odd_disease * joint_prob1 /joint_prob2

# for each patient, we first convert the hpos to the list of hpo terms such that all of them are annotated in oard

def update_hpo(hpos, method):  # input should be in "HP_..." format
    new_hpos = []
    hpo_freq = []
    for hpo in hpos:
        if hpo.replace("_", ":") in hpo_oard:
            new_hpos.append(hpo)
            hpo_freq.append(float(oard_data[oard_data['concept_code'] == hpo.replace("_", ":")]['concept_frequency'].to_string(index=False)))
        else:
            sim_score = [calc_sim(each.replace(":", "_"), hpo, method) for each in hpo_oard]
            # find the index of maximum element
            #print(sim_score)
            index_max = max(range(len(sim_score)), key=sim_score.__getitem__)
            new_hpo = hpo_oard[index_max].replace(":", "_")
            new_hpos.append(new_hpo)
            hpo_freq.append(float(oard_data[oard_data['concept_code'] == new_hpo.replace("_", ":")]['concept_frequency'].to_string(index=False)))


    return new_hpos, hpo_freq


hpo_patients = {}
freq_patients = {}
num_patients = 1


inputfile = sys.argv[1]


with open(inputfile, "r") as infile:
        for line in infile:
            hpos = line.split()[1]
            hpos = hpos.split(";")[:-1]
            update_hpos = update_hpo(hpos, 'Resnik')
            hpo_patients[num_patients] = update_hpos[0]
            freq_patients[num_patients] = update_hpos[1]
            #print(num_patients)
            num_patients+=1

hpos = []
for i in range(1,num_patients):
    hpos.extend(hpo_patients[i])

hpos = set(hpos)

disease2freq = {}
for each_hpo in hpos:
    hpo_id = domain_counts_df[domain_counts_df['concept_code'] ==each_hpo.replace("_", ":")]['concept_id'].to_string(index=False)
    params = {"dataset_id" : "2", "concept_id" :hpo_id}
    response = requests.get(url + '/frequencies/mostFrequency',params= params, verify=False)
    response_json = response.json()
    result_df = pd.DataFrame(response_json['results'])
    for each_disease in oard_disease:
        if len(result_df) ==0 or each_disease not in list(result_df['concept_id_2']):
            disease2freq[(each_hpo.replace(":", "_"), each_disease)] = -1
        else:
            disease2freq[(each_hpo.replace(":", "_"), each_disease)] = float(result_df[result_df['concept_id_2']==each_disease]['concept_frequency'].to_string(index=False))


rank_patients = {}
for i in hpo_patients.keys():
    dis_raks1 = []
    for dis in oard_disease:
        dis_raks1.append(calc_odd_oard(dis, hpo_patients[i], freq_patients[i], 0.01,disease2freq))
    rank_patients[i] = dis_raks1

    print("finished patient: "+str(i))



