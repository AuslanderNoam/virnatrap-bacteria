import sys
from scipy.stats import binom_test
from collections import defaultdict
#tcga file 1, tcga file 2, tcga file 3, gtex file, gtex index, abundance threshold

import numpy as np

def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    
def count_contigs(tokens):
    i = 1
    for token in tokens[7:]:
        if "]" in token:
            return i
        i += 1

locations = ["Gastroesophageal Junction", "Mucosa", "Muscularis"]
gtexIndex = open(sys.argv[5], "r")
gtexIndex.readline()
locDict = {}
for line in gtexIndex:
    tokens = line.split("\t")
    sample = tokens[0]
    loc = tokens[13]
    for i in range(len(locations)):
        if locations[i] in loc:
            locDict[sample] = i + 1
            break
gtexIndex.close()


gtexFile = open(sys.argv[4], "r")
gtexSamples = [defaultdict(int) for x in range(4)]
gtexDict = {}
gtexFile.readline()
for line in gtexFile:
    tokens = line.split(",")
    case = tokens[0]
    count = count_contigs(tokens)
    patient = case.split("-")[1] #Take first 4 chars of ID, which should be unique to the patient
    name = " ".join(tokens[3].split(" ")[:1]).replace("[", "").replace("]","") #genus only
    if name not in gtexDict:
        gtexDict[name] = [defaultdict(int) for x in range(4)] 
    idx = locDict[case]
    gtexDict[name][0][patient] += count
    gtexSamples[0][patient] += count
    gtexDict[name][idx][patient] += count
    gtexSamples[idx][patient] += count
gtexFile.close()

tcgaCases = defaultdict(int)
tcgaDict = {}
for i in range(1,4):
    tcgaFile = open(sys.argv[i], "r")

    tcgaFile.readline()
    for line in tcgaFile:
        tokens = line.split(",")
        case = tokens[0]
        count = count_contigs(tokens)
        name = " ".join(tokens[3].split(" ")[:1]).replace("[", "").replace("]","") #genus only
        if name not in tcgaDict:
            tcgaDict[name] = defaultdict(int)
        tcgaDict[name][case] += count
        tcgaCases[case] += count
    tcgaFile.close()

thresh = float(sys.argv[6])
    
for count_thresh in range(1, 11):
    print("Threshold:", count_thresh)
    
    #print(",".join(("Genus", "TCGA Abundance", "GTEx Overall Abundance", "Pvalue", "Corrected Pvalue", "TCGA Greater", "GTEx Junction Abundance", "GTEx Mucosa Abundance", "GTEx Muscularis Abundance")))
    allSpecies = set(tcgaDict.keys()) | set(gtexDict.keys())
    
    #Remove below-threshold (case,genera) pairs
    allDicts = gtexSamples + [tcgaCases] + [d for sp in gtexDict.keys() for d in gtexDict[sp]] + [tcgaDict[sp] for sp in tcgaDict.keys()]
    for d in allDicts:
        items = list(d.items())
        for case, count in items:
            if count < count_thresh:
                del d[case]         

    tcgaCount = len(tcgaCases)
    gtexCounts = [len(x) for x in gtexSamples]
    print("Overall TCGA and GTEx counts:")
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
            dir = 0
        else:
            p = binom_test(tcgaPos, n=tcgaCount, p=adjGtexP, alternative='greater') #anchored on gtex p
            dir = 1
        pValues.append(p)
        linesToPrint.append([sp, tcgaFrac, gtexFracs[0], p, -1, dir])

    pAdjust = p_adjust_bh(pValues)
    count_esca = 0
    count_gtex = 0
    count_non = 0
    for i in range(len(linesToPrint)):
        tp = linesToPrint[i]
        tp[4] = pAdjust[i]
        if tp[4] < 0.05:
            if tp[5]:
                count_esca += 1
            else:
                count_gtex += 1
        else:
            count_non += 1

    print("Signficant + in ESCA:", count_esca)
    print("Signficant - in ESCA:", count_gtex)
    print("Not significant:     ", count_non)
    print("\n")






