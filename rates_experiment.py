import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import re
from collections import Counter
from ete3 import Tree
from tabulate import tabulate

from lingdata import database
from lingdata.categorical import CategoricalData

import code.pythia as pythia
import code.mptp as mptp
import code.raxmlng as raxmlng
import code.util as util



def run_raxml_ng(df):
    for (i, row) in df.iterrows():
        raxmlng.run_inference_adaptive(row["msa_paths"]["bin"], "BIN+G", util.prefix(results_dir, row, "gamma_adaptive", "bin"))


def evaluate(df):
    average_group_sizes = [[] for _ in range(4)]
    pos_array = [[0 for __ in range(4)] for _ in range(4)]
    analyzed = 0
    for i, row in df.iterrows():
        r = []
        averages = []
        prefix = util.prefix(results_dir, row, "gamma_adaptive", "bin")
        alpha = raxmlng.alpha(prefix)
        if alpha >= 20: #low heterogeneity
            continue
        srlhs = raxmlng.siterate_lhs(prefix)
        weights_rates = raxmlng.weights_rates(prefix)
        if len(srlhs) == 0:
            continue
        #print(row["ds_id"])
        #print("alpha: " + str(alpha))
        best_rates = []
        for rlhs in srlhs:
            best_llh = max(rlhs)
            best_indices = []
            for r_i in range(4):
                if rlhs[r_i] == best_llh:
                    best_indices.append(r_i)
            best_rates.append(best_indices)
        buckets = [[] for _ in range(4)]
        site_group_sizes = row["site_group_sizes"]
        assert(len(site_group_sizes))
        for i, group_size in enumerate(site_group_sizes):
            for r_i in best_rates[i]:
                buckets[r_i].append(group_size)
        for i, bucket in enumerate(buckets):
            if len(bucket) == 0:
                a = 0
            else:
                a = sum(bucket) / len(bucket)
            r.append([i, weights_rates[i][1], len(bucket), a])
            averages.append(a)
            average_group_sizes[i].append(a)
        #print(tabulate(r, tablefmt="pipe", floatfmt=".2f", headers = ["rate_index", "rate", "number of sites", "average site group size"]))
        L = [(averages[r_i],r_i) for r_i in range(4)]
        L.sort()
        _,permutation = zip(*L)
        for r_i in range(4):
            pos_array[r_i][permutation.index(r_i)] += 1
        analyzed +=1
    #for rate_i in range(4):
    #    plt.hist(average_group_sizes[rate_i], alpha=0.5, label=str(rate_i))
    #    plt.xlim(0, 30)
    #    plt.ylim(0, 30)
    #    plt.legend(loc='upper right')
    #    plt.xlabel("Average group size")
    #    plt.ylabel("Number of datasets")
    #    plt.savefig(os.path.join(plots_dir, "site_groups_rates_" + str(rate_i) + ".png"))
    #    plt.clf()
    fig,ax = plt.subplots()
    y_old = [0 for r_i in range(4)]
    x = range(4)
    for r_i in range(4):
        ax.bar(range(4), pos_array[r_i], bottom=y_old, label = str(r_i))
        y_old = [y_old[r_j] + pos_array[r_i][r_j] for r_j in range(4)]
    plt.legend()
    #plt.xlabel("Position in ordering by average site group size")
    #plt.ylabel("Number of datasets")
    plt.savefig(os.path.join(plots_dir, "high_het_pos.png"))
    plt.clf()


raxmlng.adaptive_exe_path = "./bin/raxml-ng-adaptive"
mptp.exe_path = "./bin/mptp"
pythia.raxmlng_path = "./bin/raxml-ng"
pythia.predictor_path = "predictors/latest.pckl"
config_path = "lingdata_modeltesting_familyfull_config.json"

database.read_config(config_path)
#database.compile()
df = database.data()
results_dir = os.path.join("data/results", "familyfull")
plots_dir = os.path.join("data/results", "plots")
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
#run_raxml_ng(df)
evaluate(df)
