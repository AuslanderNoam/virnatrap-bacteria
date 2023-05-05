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

DESCRIPTION="ToDo"


def blastx_contig(args):
    file = args[0]
    name = args[1]
    nm = os.getcwd()+"/"+name+"/" + \
         file.replace('_gdc_realn_rehead_contigs.txt', '.txt').split('/')[-1]
    subprocess.run('blastx -db $BLASTDB/bacprot0.fa -query ' + file + ' -out ' + nm + ' -outfmt 6 -evalue 0.00001',
                   shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--input_dir_contigs", type=str, help="input directory", required=True
    )
    parser.add_argument(
        "--name", type=str, help="name of project", required=False
    )
    args = parser.parse_args()

    input_dir_contigs = args.input_dir_contigs
    name = args.name

    fls = glob.glob(input_dir_contigs + '*.txt')
    try:
        fls.remove("/wistar/dbgap-mutporg/download_tools/GTEx/ESOPHAGUS/Results/assembly_output2/GTEX-OIZF-0626-SM-2I5GT_contigs.txt") #This is the in progress ESO sample
    except ValueError:
        pass     
    pool = mp.Pool(processes=48)
    os.mkdir(os.getcwd() + "/" + name)
    pool.map(blastx_contig, [[f,name] for f in fls])