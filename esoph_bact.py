import sys
from scipy.stats import binom_test
#tcga file 1, tcga file 2, tcga file 3, gtex file, gtex index, threshold

import numpy as np

def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

locations = ["Gastroesophageal Junction", "Mucosa", "Muscularis"]
gtexIndex = open(sys.argv[5], "r")
gtexIndex.readline()
locDict = {}
for line in gtexIndex:
    tokens = line.split("\t")
    sample = tokens[0]
    loc = tokens[13]
    for i in range(len(locations)):
        if locations[i]  in loc:
            locDict[sample] = i + 1
            break
gtexIndex.close()

gtexFile = open(sys.argv[4], "r")
gtexSamples = [set(), set(), set(), set()]
gtexDict = {}
gtexFile.readline()
for line in gtexFile:
    tokens = line.split(",")
    case = tokens[0]
    patient = case.split("-")[1] #Take first 4 chars of ID, which should be unique to the patient
    name = " ".join(tokens[3].split(" ")[:1]).replace("[", "").replace("]","") #genus only - modify [:1] for full names
    if name not in gtexDict:
        gtexDict[name] = [set(), set(), set(), set()]   
    idx = locDict[case]
    gtexDict[name][0].add(patient)
    gtexSamples[0].add(patient)    
    gtexDict[name][idx].add(patient)
    gtexSamples[idx].add(patient)
gtexFile.close()

tcgaCases = set()
tcgaDict = {}
for i in range(1,4):
    tcgaFile = open(sys.argv[i], "r")

    tcgaFile.readline()
    for line in tcgaFile:
        tokens = line.split(",")
        case = tokens[0]
        name = " ".join(tokens[3].split(" ")[:1]).replace("[", "").replace("]","") #genus only - modify [:1] for full names
        if name not in tcgaDict:
            tcgaDict[name] = set()
        tcgaDict[name].add(case)
        tcgaCases.add(case)
    tcgaFile.close()


thresh = float(sys.argv[6])
print(",".join(("Genus", "TCGA Abundance", "GTEx Overall Abundance", "Pvalue", "Corrected Pvalue", "TCGA Greater", "GTEx Junction Abundance", "GTEx Mucosa Abundance", "GTEx Muscularis Abundance")))
allSpecies = set(tcgaDict.keys()) | set(gtexDict.keys())
tcgaCount = len(tcgaCases)
gtexCounts = [len(x) for x in gtexSamples]
print(tcgaCount, gtexCounts, file=sys.stderr)
linesToPrint = []
pValues = []
for sp in allSpecies:
    tcgaPos = len(tcgaDict[sp]) if sp in tcgaDict else 0
    gtexPos = len(gtexDict[sp][0]) if sp in gtexDict else 0 
    tcgaFrac = tcgaPos/tcgaCount
    gtexFracs = [len(gtexDict[sp][i])/gtexCounts[i] for i in range(len(gtexCounts))] if sp in gtexDict else [0 for i in range(len(gtexCounts))]
    if max(tcgaFrac, gtexFracs[0]) < thresh:
        continue
    adjGtexP = gtexFracs[0]
    #p = 0 or 1 exactly gives a p-val of 0 no matter what from the binom test
    adjGtexP = max(0.0001, adjGtexP)
    adjGtexP = min(0.9999, adjGtexP)
    if tcgaFrac < gtexFracs[0]:
        p = binom_test(tcgaPos, n=tcgaCount, p=adjGtexP, alternative='less')
        dir = '0'
    else:
        p = binom_test(tcgaPos, n=tcgaCount, p=adjGtexP, alternative='greater') #anchored on gtex p
        dir = '1'
    pValues.append(p)
    linesToPrint.append([sp, str(tcgaFrac), str(gtexFracs[0])] + [str(p), '-1', dir] + [str(x) for x in gtexFracs[1:]])
        
pAdjust = p_adjust_bh(pValues)
for i in range(len(linesToPrint)):
    tp = linesToPrint[i]
    tp[4] = str(pAdjust[i])
    print(",".join(tp))





