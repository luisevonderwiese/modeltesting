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
        all_alphas.append(raxmlng.get_partitioning_alphas(util.prefix(results_dir, row, "partitioning", "bin_BIN+G_x")))
    max_part = max([max(alphas.keys()) for alphas in all_alphas if alphas != {}])
    alphas_matrix = []
    labels = ["[0,20]", "[20,90]", "[90, 100]"]
    for idx in range(max_part):
        idx_alphas = [0, 0, 0]
        for alphas in all_alphas:
            if idx in alphas:
                alpha = alphas[idx]
                if alpha < 20:
                    idx_alphas[0] += 1
                elif alpha < 90:
                    idx_alphas[1] += 1
                else:
                    idx_alphas[2]+= 1
        alphas_matrix.append(idx_alphas)
    fig,ax = plt.subplots(figsize=(15, 10))
    x = range(len(alphas_matrix))
    y_old = [0 for el in x]
    for num in range(3):
        y_new = []
        for idx_alphas in alphas_matrix:
            y_new.append(idx_alphas[num])
        ax.bar(x, y_new, bottom=y_old, label = labels[num])
        for i in x:
            y_old[i] = y_old[i] + y_new[i]
    ax.legend()
    plt.xlabel("Site Group Size")
    plt.ylabel("Number of Datasets")
    plt.savefig(os.path.join(plots_dir, "stacked_alpha_ranges.png"))
    plt.clf()

raxmlng.exe_path = "./bin/raxml-ng"
config_path = "lingdata_modeltesting_familyfull_partitioning.json"
database.read_config(config_path)
#database.compile()
df = database.data()
results_dir = os.path.join("data/results/partitioning")
plots_dir = os.path.join("data/results/partitioning_plots")
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)

#run_raxml_ng(df)
analyze_alphas(df)
