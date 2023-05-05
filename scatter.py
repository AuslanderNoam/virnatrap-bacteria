import sys
import matplotlib.pyplot as plt
#CSV, cancer list, normal list

def parseSpecies(file):
    species = set()
    inSpecies = open(file)
    for line in inSpecies:
        sp = line.strip()
        if len(sp) > 0:
            species.add(sp)
    inSpecies.close()
    return species

cancerSpecies = parseSpecies(sys.argv[2])#.intersection(allSpecies)
healthySpecies = parseSpecies(sys.argv[3])#.intersection(allSpecies)

cols = []
tcga = []
gtex = []
inFile = open(sys.argv[1])
inFile.readline()
for line in inFile:
    tokens = line.split(",")
    idx = 0
    col = 'gray'
    n = tokens[0]
    t = float(tokens[1])
    g = float(tokens[2])
    if abs(t-g) > 0.5:
        print(n, g, t)
    if n in cancerSpecies:
        col = 'red'
        idx = len(tcga)
    elif n in healthySpecies:
        col = 'blue'
        idx = len(tcga)
    
    tcga.insert(idx, t)
    gtex.insert(idx, g)
    cols.insert(idx, col)
inFile.close()

plt.figure(figsize=(6,6))
plt.scatter(gtex, tcga, c=cols, alpha=0.5)
plt.xlabel("GTEx Abundance", size='large')
plt.ylabel("TCGA Abundance", size='large')
plt.savefig("genus_scatter.svg")
