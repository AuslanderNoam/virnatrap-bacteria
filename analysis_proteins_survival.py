import pandas as pd
import numpy as np
import seaborn as sns
import scipy
import lifelines
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import ranksums
from scipy import stats
plt.rcParams["font.family"] = "Arial"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.facecolor'] = 'none'
findi = lambda x,v:[i for i in range(len(x)) if x[i]==v]

def plot_km(BAC,OST,OSS, mkplt = True, THR=0):
    plt.figure()
    if BAC.name is None:
        BAC.name = 'prots'
    death = [i=='1:DECEASED' or i=='1:DEAD WITH TUMOR' for i in OSS]
    dat = pd.DataFrame({'surv': OST, 'death': death, 'val': BAC})
    lr = lifelines.statistics.logrank_test(dat[dat.val <= THR].surv, dat[dat.val > THR].surv,
                                           event_observed_A=dat[dat.val <= THR].death,
                                           event_observed_B=dat[dat.val > THR].death)
    kmf = KaplanMeierFitter()
    kmf.fit(dat[dat.val <= THR].surv, event_observed=dat[dat.val <= THR].death, label=BAC.name + " negative")
    s1 = kmf.median_survival_time_
    if mkplt:
        ax=kmf.plot_survival_function()

    kmf.fit(dat[dat.val > THR].surv, event_observed=dat[dat.val > THR].death, label=BAC.name +" positive")
    s2 = kmf.median_survival_time_
    if mkplt:
        kmf.plot_survival_function(ax=ax)
        ax.set_ylabel('survival probability')
        ax.set_xlabel('time')
        mytext = "P=%.2e" % (lr.p_value)
        ax.text(0.1, 0.1, mytext)
    return lr.p_value,s1/s2

##Load ESCA and ESO protein mapping results
esca = pd.read_csv('bac_table_esca.csv', index_col = 'case_id')
eso = pd.read_csv('bac_table_eso.csv', index_col = 'case_id')

#read and merge clinical data
clin1 = pd.read_csv('esca_tcga_pan_can_atlas_2018/data_clinical_sample.csv')
clin1.index = clin1['PATIENT_ID']

clin2 = pd.read_csv('esca_tcga_pan_can_atlas_2018/data_clinical_patient.csv')
clin2.index = clin2['PATIENT_ID']

clin = clin1.join(clin2,lsuffix='_p',)

#get samples in both tables and orginize similarly
inters = list(set(clin.index)&set(esca.index))

#get orginized data
esca = esca.loc[inters]
clin = clin.loc[inters]

##get the full names of these proteins from the fasta database
ff = open('bacprot0.fa')
ll=ff.readlines()
nn = [i for i in ll if '>' in i]
dd = {i.split(' ')[0].replace('>',''):''.join(i.split(' ')[1:]) for i in nn}

#select bacterial proteins to be tested for survival - proteins with at least 5 sampeles with/without the proteins
ll=list(esca.keys())
ss=esca.sum()
bacn = [i for i in ll if ss[i]>=5 and ss[i]<=(esca.shape[0]-5)]

OST = clin['OS_MONTHS']
OSS = clin['OS_STATUS']

#get individual bacteria associated with OS through log-rank pvalue
pvl = []
hr=[]
for b in bacn:
    BAC = esca[b]
    p, h = plot_km(BAC, OST, OSS, mkplt=False)
    hr.append(h)
    pvl.append(p)

##FDR-correct for multiple hypotheses
rejected_ps,pvlc_os = fdrcorrection(pvl, alpha=0.05)

##29 significant hits: ALL associated with poor survival (hr>0)
nms_poor0 = [bacn[i] for i in range(len(pvl)) if pvlc_os[i]<0.05 and hr[i]>1]

OST = clin['DSS_MONTHS']
OSS = clin['DSS_STATUS']
#get individual bacteria associated with OS through log-rank pvalue
pvl2 = []
hr2=[]
for b in bacn:
    BAC = esca[b]
    p,h=plot_km(BAC, OST, OSS, mkplt=False)
    hr2.append(h)
    pvl2.append(p)

##FDR-correct for multiple hypotheses
rejected_dss,pvlc_dss = fdrcorrection(pvl2, alpha=0.05)

##29 significant hits: 24 associated with poor survival (hr>0)
nms_poor = [bacn[i] for i in range(len(pvl2)) if pvlc_dss[i]<0.05 and hr2[i]>1]
##29 significant hits: 5 associated with better survival (hr<0)
nms_better = [bacn[i] for i in range(len(pvl2)) if pvlc_dss[i]<0.05 and hr2[i]<1]

##get all survival results and save
rest = pd.DataFrame({'microbial_protein_name':bacn,'os_ratio':hr,'os_q':pvlc_os,'dss_ratio':hr2,'dss_q':pvlc_dss}).to_csv('Results/survival_analysis.csv')

##plot these five KMs (non-redundant, interesting, significant through OS and DSS)
prots_km = list(set(nms_poor0)&set(nms_poor))
prots_km = [prots_km[i] for i in [2,4,5,6,7]]

for i in range(len(prots_km)):
    prt = esca[prots_km[i]]
    plot_km(prt, OST, OSS, mkplt=True)


##survival analysis considering any Fe protein
bp = esca[prt]
z = bp.sum(1)
plot_km(z, OST, OSS, mkplt=True)

exp = pd.read_csv('esca_tcga_pan_can_atlas_2018/tcga_exp.csv',index_col='Hugo_Symbol')
exp=exp[inters]
prt=['WP_006680945.1','WP_002532908.1','WP_131625607.1']
pve=[];pve2=[];m1=[];m2=[]
for i in range(len(exp)):
    ge = exp.iloc[i]
    r, pv = ranksums(ge[z > 0], ge[z == 0],alternative='greater')
    pve.append(pv)
    r, pv = ranksums(ge[z > 0], ge[z == 0], alternative='less')
    pve2.append(pv)
    m1.append(np.median(ge[z > 0]))
    m2.append(np.median(ge[z == 0]))

pd.DataFrame([exp.index[i] for i in range(len(exp.index)) if pve[i]<0.05 and m1[i]>0.2 and m2[i]<0]).to_csv('Results/genes_upregulatde.csv')
pd.DataFrame([exp.index[i] for i in range(len(exp.index)) if pve2[i]<0.05 and m1[i]<-0.2 and m2[i]>0]).to_csv('Results/genes_downregulatde.csv')

##Iron bacterial proteins are associated with up regulation of ferrosig,infections, oxphos,endocytosis
##only ferrosig genes are also associated with poor survival
# ferroptosis upregulated genes
ferrosig=['MAP1LC3B','MAP1LC3B2','VDAC2','SAT2','SAT1','FTL']
# infection upregulated genes
infection = ['ATP6V1H', 'ATP6V1G1', 'DCTN3', 'MAP2K1', 'DYNLL1', 'TUBB2A', 'NDUFB11', 'DYNLT1', 'COX6B1', 'ATP6V0D1', 'NDUFB2', 'HBEGF', 'NDUFB9', 'NDUFB3', 'COX4I1', 'NDUFB5', 'ATP6V0C', 'NDUFB8', 'UQCRC1', 'COX7A2L', 'CSK', 'GAPDH', 'ACTB', 'NDUFS2', 'COX6A1', 'NDUFC2', 'NDUFA3', 'COX8A', 'NDUFC1', 'NDUFS3', 'NDUFB7', 'BAX', 'VPS16', 'ARPC3', 'ATP6V1E1', 'NDUFA1', 'HRAS', 'ATP6V0B', 'NDUFA13', 'TRADD', 'ACTR3', 'ARPC5L', 'PPA1', 'COX5A', 'ARL8A', 'COX17', 'ARF1', 'KLC2', 'NDUFS6', 'COX7B', 'KPNA3', 'ATP6V1F', 'ATP6AP1', 'IL6', 'IKBKG']
# oxphos upregulated genes
oxphos=['NDUFB9', 'NDUFA13', 'COX7B', 'NDUFB8', 'NDUFB7', 'NDUFB11', 'NDUFB5', 'COX4I1', 'COX17', 'NDUFB3', 'NDUFB2', 'COX6A1', 'COX5A', 'ATP6V1H', 'ATP6V1E1', 'ATP6V1F', 'COX8A', 'ATP6V1G1', 'ATP6V0B', 'ATP6AP1', 'NDUFA3', 'NDUFA1', 'NDUFC2', 'NDUFC1', 'COX6B1', 'COX7A2L', 'PPA1', 'NDUFS6', 'UQCRC1', 'NDUFS3', 'NDUFS2', 'ATP6V0D1', 'ATP6V0C']
# endocytosis upregulated genes
endocytosis=['ACTR3', 'SH3GLB2', 'ARF1', 'VPS4A', 'ARPC5L', 'VPS26A', 'AP2A1', 'SNF8', 'EPS15L1', 'EPN1', 'HGS', 'ARPC3', 'CHMP1A', 'AP2S1', 'CHMP2A', 'CHMP4A', 'HRAS', 'RAB8A', 'AP2M1', 'SPG21']

patd = pd.DataFrame({'bacprot':z>0,'ferroptosis':exp.loc[ferrosig].mean(),'bacterial_infection':exp.loc[infection].mean(),'oxphos':exp.loc[oxphos].mean(),'endocytosis':exp.loc[endocytosis].mean()})

# Make barplots for Figure 4b
ax1 = plt.subplot(1, 4, 1)
sns.boxplot(data=patd, x="bacprot", y="ferroptosis", hue="bacprot",ax=ax1)
ax2 = plt.subplot(1, 4, 2)
sns.boxplot(data=patd, x="bacprot", y="bacterial_infection", hue="bacprot",ax=ax2)
ax3 = plt.subplot(1, 4, 3)
sns.boxplot(data=patd, x="bacprot", y="oxphos", hue="bacprot",ax=ax3)
ax4 = plt.subplot(1, 4, 4)
sns.boxplot(data=patd, x="bacprot", y="endocytosis", hue="bacprot",ax=ax4)

##survival analysis considering Fe signature
val = stats.zscore(exp.loc[ferrosig], axis=1).mean()
plot_km(val, OST, OSS, mkplt=True)

# Make heatmaps in Figure 4a
gg = list(set(ferrosig+infection+oxphos+endocytosis))
xx = exp.loc[gg]
x1=xx[z[z>0].index]
x0=xx[z[z==0].index]
xx2=xx[list(z[z>0].index)+list(z[z==0].index)]

row_colors=['r' for i in range(x1.shape[1])]+['b' for i in range(x0.shape[1])]
col_colors=['silver' for i in gg]
for i in range(len(gg)):
    if gg[i] in ferrosig:
        col_colors[i] = 'silver'
    if gg[i] in oxphos and gg[i] in infection:
        col_colors[i] = 'skyblue'
    if gg[i] in oxphos and gg[i] not in infection:
        col_colors[i] = 'steelblue'
    if gg[i] in endocytosis and gg[i] not in infection:
        col_colors[i] = 'lightcoral'
    if gg[i] in endocytosis and gg[i]  in infection:
        col_colors[i] = 'blueviolet'
    if gg[i] in infection and gg[i] not in endocytosis and gg[i] not in oxphos:
        col_colors[i] = 'orange'

#make heatmap
sns.clustermap(xx2.T,row_cluster=False,xticklabels=1,cmap="PiYG",vmin=-3, vmax=3,row_colors=row_colors,col_colors=col_colors)# prot_high_sig_os=['WP_002477833.1', 'WP_002477847.1', 'WP_024181148.1','WP_033636612.1','WP_046466356.1','WP_131625607.1']
