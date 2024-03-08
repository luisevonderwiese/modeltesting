import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import re
from collections import Counter
from ete3 import Tree

from lingdata import database
from lingdata.categorical import CategoricalData

import code.pythia as pythia
import code.mptp as mptp
import code.raxmlng as raxmlng
import code.util as util


def run_raxml_ng(df):
    for (i, row) in df.iterrows():
        raxmlng.run_inference(row["msa_paths"]["bin"], row["partition_paths"]["bin_BIN+G_x"], util.prefix(results_dir, row, "partitioning", "bin_BIN+G_x"))

def analyze_alphas(df):
    all_alphas = []
    for (i, row) in df.iterrows():
        prefix = util.prefix(results_dir, row, "partitioning", "bin_BIN+G_x")
        all_alphas.append(raxmlng.get_partitioning_alphas(prefix))
    max_part = max([max(alphas.keys()) for alphas in all_alphas if alphas != {}])
    alphas_matrix = [[] for _ in range(min(26, max_part))]
    for idx in range(min(26, max_part)): #cutoff because there are only few large site groups in the data
        for alphas in all_alphas:
            if idx in alphas:
                alphas_matrix[idx].append(alphas[idx])
    nans = [float('nan'), float('nan')] # requires at least 2 nans
    fig,ax = plt.subplots(figsize=(30, 10))
    ax.violinplot([alphas or nans for alphas in alphas_matrix])
    plt.xlabel("Site Group Size")
    plt.ylabel("Alpha")
    plt.savefig(os.path.join(plots_dir, mode + "_alpha_violinplot.png"))
    plt.clf()

mode = "familysplit"
raxmlng.exe_path = "./bin/raxml-ng"
config_path = "lingdata_modeltesting_" + mode + "_partitioning.json"
database.read_config(config_path)
#database.compile()
df = database.data()
results_dir = os.path.join("data/results/" + mode + "_partitioning")
plots_dir = os.path.join("data/results/partitioning_plots")
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)

#run_raxml_ng(df)
analyze_alphas(df)
