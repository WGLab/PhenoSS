#!/usr/bin/python

from __future__ import print_function # load print function in python3
from collections import defaultdict
import numpy as np

import ssmpy, sys

ssmpy.ssm.mica = True
ssmpy.ssm.intrinsic = True

ssmpy.semantic_base("hp.db")

# set up auto dictionary function
def auto_dict():
    return defaultdict(auto_dict)

inputfile = sys.argv[1]

patient = sys.argv[2]



dis1 = auto_dict()
dis2 = auto_dict()
comparedone = auto_dict()

with open(inputfile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        tmpinf = line.split("\t")
        dis1[tmpinf[0]] = tmpinf[1]


with open(inputfile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        tmpinf = line.split("\t")
        dis2[tmpinf[0]] = tmpinf[1]

for pat1 in dis1:
    for pat2 in dis2:
        comparedone[pat1][pat2] = 0
        comparedone[pat2][pat1] = 0

outFile = patient + "_sim"
OUT = open(outFile, 'w')
        
for pat1 in dis1:
    if pat1 != patient:
        continue
    for pat2 in dis2:
    
        hpo1 = dis1[pat1]
        hpo2 = dis2[pat2]

        #print(hpo1+"\t"+hpo2)
        
        hpo1 = hpo1.split(";")
        hpo2 = hpo2.split(";")

        totaldiff = 0
        numcompare = 0
        
        for hp1 in hpo1:
            #print(pat1 + "\t" +hp1)
            alldiff = []
            numcompare += 1
            for hp2 in hpo2:
                e1 = ssmpy.get_id(hp1)
                e2 = ssmpy.get_id(hp2)
                diff_resnik = ssmpy.ssm_resnik(e1,e2)
                alldiff.append(diff_resnik)
                out = pat1+"\t"+pat2+"\t"+ hp1 + "\t" + hp2 + "\t" + str(diff_resnik)
                #print(out)
            maxdiff = np.amax(alldiff)
            totaldiff += maxdiff
                
        if numcompare == 0:
            out = pat1 + "\t" + pat2 + "\tNA"
        else:
            diff = totaldiff/numcompare
            out = pat1 + "\t" + pat2 + "\t" + str(diff)
            #comparedone[pat1][pat2] = diff
            #comparedone[pat2][pat1] = diff
        print(pat2)
        OUT.write(out+"\n")
OUT.close()
