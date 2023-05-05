import pandas as pd
import numpy as np
import seaborn as sns
import scipy
import lifelines
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import ranksums

plt.rcParams["font.family"] = "Arial"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.facecolor'] = 'none'

findi = lambda x,v:[i for i in range(len(x)) if x[i]==v]


##get the full names of these proteins from the fasta database
ff = open('bacprot0.fa')
ll=ff.readlines()
nn = [i for i in ll if '>' in i]
dd = {i.split(' ')[0].replace('>',''):''.join(i.split(' ')[1:]) for i in nn}


##bacterial proteins
esca = pd.read_csv('bac_table_esca.csv', index_col = 'case_id')
eso = pd.read_csv('bac_table_eso.csv', index_col = 'case_id')

##Combine samples from the same GTEx patient
eso['patient_id'] = [x.split("-")[1] for x in eso.index]
eso = eso.groupby(eso['patient_id']).aggregate(max)

##get proteins found in either and fill with zeros
all_prots = set(list(esca.columns)+list(eso.columns))
ab1 = [i for i in all_prots if i not in esca.columns]
ab2 = [i for i in all_prots if i not in eso.columns]

for i in range(len(ab1)):
    esca[ab1[i]] = [0 for i in range(esca.shape[0])]

for i in range(len(ab2)):
    eso[ab2[i]] = [0 for i in range(eso.shape[0])]
esca = esca[eso.columns]

#get frequencies of proteins
s1 = esca.sum()
freq_esca = [s1[i]/esca.shape[0] for i in range(esca.shape[1])]

s2 = eso.sum()
freq_eso = [(s2[i])/eso.shape[0] for i in range(eso.shape[1])]

freq_ratios = [(freq_esca[i]+1)/(freq_eso[i]+1) for i in range(len(freq_eso))]
freq_diffs = [freq_esca[i]-freq_eso[i] for i in range(len(freq_eso))]
dfreq = {eso.columns[i]:[freq_esca[i],freq_eso[i]] for i in range(len(freq_esca))}

#save data of protein frequency
results_freqs = pd.DataFrame({'microbial_protein_name':list(esca.columns),'ESCA_freq':freq_esca,'ESO_freq':freq_eso,'freq_score':freq_ratios,'freq_diff':freq_diffs})
results_freqs = results_freqs.sort_values(by=['freq_diff'])
results_freqs.to_csv('Results/protein_frequencies_fix.csv')

#proteins to plot - picked from top/bottom rations
include_1 = ['WP_011114055.1','WP_073668305.1','WP_087674308.1', 'WP_185964076.1','WP_011114056.1'] ##manually take unique characterized ones higher in GTEx
include_2 = ['WP_001359455.1','WP_023147606.1','WP_134326155.1','WP_202116303.1','WP_087694260.1','WP_134327566.1','WP_000063280.1','WP_000453580.1'] ##manually take unique characterized ones higher in TCGA
include_prots = include_1+include_2


##load bacterial species
esca_b = pd.read_csv('combined_ESCA.csv', index_col = 'case_id')
eso_b = pd.read_csv('combined_ESO.csv', index_col = 'case_id')

##proteins with umnapped species
esca_b['bacteria_name']= [i.split(' ')[0] if not isinstance(i, float) else 'other' for i in list(esca_b.bacteria_name)]
eso_b['bacteria_name']= [i.split(' ')[0] if not isinstance(i, float) else 'other' for i in list(eso_b.bacteria_name)]

#map selected proteins to bacteria sepcies
bacsc=[]
bacsn=[]
for i in range(len(include_prots)):
    bacsc = bacsc + list(esca_b[esca_b.bac == include_prots[i]].bacteria_name)
    bacsn = bacsn + list(eso_b[eso_b.bac == include_prots[i]].bacteria_name)

#select species to plot - at least 5% of sample size
v1 = [i[0] for i in pd.DataFrame(bacsc).value_counts()[:14].index] ##manual at least 5% of samples
v2 = [i[0] for i in pd.DataFrame(bacsn).value_counts()[:19].index]##manual at least 5% of samples


# get the read,genus relations for normal and cancer, to plot F3A
genus = list(set(v1+v2))
c1=[];c2=[]
for i in range(len(include_prots)):
    z = esca_b[esca_b.bac == include_prots[i]]
    v = eso_b[eso_b.bac == include_prots[i]]
    cc1=[]
    cc2=[]
    for j in range(len(genus)):
        cc1.append(len(z[z.bacteria_name==genus[j]]))
        cc2.append(len(v[v.bacteria_name == genus[j]]))
    c1.append(cc1) #cancer
    c2.append(cc2) #normal

c1 = np.array(c1)
c2 = np.array(c2)

# reorder data
g = sns.clustermap(c1)
g.dendrogram_col.reordered_ind
col = [genus[i] for i in g.dendrogram_col.reordered_ind]
row = [include_prots[i] for i in g.dendrogram_row.reordered_ind]


# save to make source data of Fig3 panel A
x = pd.DataFrame(c1)
x.index = include_prots
x.columns = genus
x = x.loc[row]
x = x[col]
x.to_csv('Results/canc_table_fix.csv')

x = pd.DataFrame(c2)
x.index = include_prots
x.columns = genus
x = x.loc[row]
x = x[col]
x.to_csv('Results/norm_table_fix.csv')


# bar plots of the frequencies
freqs = pd.DataFrame({'ESCA_freq':[dfreq[i][0] for i in x.index],'ESO_freq':[dfreq[i][1] for i in x.index]},index=[x.index])

freqs['ESCA_freq'].plot.bar()
freqs['ESO_freq'].plot.bar()
