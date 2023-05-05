import os
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
import json
import re
import argparse
import random
import subprocess
import multiprocessing as mp
from multiprocessing import freeze_support
from os.path import exists
import ast


DESCRIPTION="ToDo"

def extract_combine(args):
    path_csv = args[0]
    name = args[1]
    blastoutpath = os.getcwd() + "/" + name + "/"
    files = glob.glob(blastoutpath+'*.txt')
    col_names = ["id", "bac", "score", "coverage", "blst1", "blst2", "blst3", "blst4", "blst5", "blst6", "blst7", "blst8"]

    ##1. read all the blastx outputs and combine to one table
    flag = 1
    for file in files:
        x = pd.read_csv(file, sep='\t', header=None, names=col_names)
        x2 = x[x.score > 95]
        id1 = list(x2['id'])
        x2['id'] = [i.split('_bact:')[0] for i in id1]
        x2['sample'] = file.split('/')[-1].split('.')[0].split('_')[0]
        if flag:
            combined = x2
            flag = 0
        else:
            combined = pd.concat([combined, x2])

    combined.to_csv('combined.csv')
    ##2. add sample and contig information for every contig mapped with blastx
    resc = pd.read_csv(path_csv)
    resc.index = resc['sample']

    lc = list(resc['contigs'])
    dictc = {}
    for i in range(len(lc)):
        cc = ast.literal_eval(lc[i])
        for co in cc:
            dictc[co] = i

    smp = list(combined['sample'])
    contigs = list(combined['id'])

    #terrible pandas nonsense
    combined['project_id'] = ['NA' for i in smp]
    combined['case_id'] = ['NA' for i in smp]
    combined['sample_type'] = ['NA' for i in smp]
    combined['bacteria_name'] = ['NA' for i in smp]
    combined['bacteria_acc'] = ['NA' for i in smp]
    #trusting that the intersect of samples is complete otherwise NAN --- VERIFY
    for i in range(len(smp)):
        if smp[i] in resc.index:
            lc = resc.loc[smp[i]].iloc[0]
            combined['project_id'].iloc[i] = lc.project_id
            combined['case_id'].iloc[i] = lc.case_id
            combined['sample_type'].iloc[i] = lc.sample_type
            if contigs[i] in dictc.keys():
                combined['bacteria_name'].iloc[i] = resc.iloc[dictc[contigs[i]]]['seq_name']
                combined['bacteria_acc'].iloc[i] = resc.iloc[dictc[contigs[i]]]['accession']

    combined.to_csv('combined_'+name+'.csv')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--path_csv", type=str, help="path_csv", required=True
    )
    parser.add_argument(
        "--name", type=str, help="name of project", required=False
    )
    args = parser.parse_args()

    path_csv = args.path_csv
    name = args.name


    extract_combine([path_csv,name])

