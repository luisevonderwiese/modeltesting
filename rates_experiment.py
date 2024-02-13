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
        raxmlng.run_inference_adaptive(row["msa_paths"]["bin"], "BIN+G", util.prefix(results_dir, row, "gamma_adaptive", "bin"))


def evaluate(df):
    average_group_sizes = [[] for _ in range(4)]
    analyzed = 0
    for i, row in df.iterrows():
        prefix = util.prefix(results_dir, row, "gamma_adaptive", "bin")
        srlhs = raxmlng.siterate_lhs(prefix)
        if len(srlhs) == 0:
            continue
        best_rates = []
        for rlhs in srlhs:
            best_rates.append(rlhs.index(max(rlhs)))
        buckets = [[] for _ in range(4)]
        site_group_sizes = row["site_group_sizes"]
        assert(len(site_group_sizes))
        for i, group_size in enumerate(site_group_sizes):
            buckets[best_rates[i]].append(group_size)
        for i, bucket in enumerate(bucktes):
            average_group_sizes[i].append(sum(bucket) / len(bucket))
        analyzed +=1
    for rate_i in range(4):
        plt.hist(average_group_sizes[rate_i], alpha=0.5, label=str(rate_i))

    plt.legend(loc='upper right')
    plt.xlabel("Average group size")
    plt.ylabel("Number of datasets")
    plt.savefig(os.path.join(plots_dir, "site_groups_rates.png"))
    plt.clf()


raxmlng.adaptive_exe_path = "./bin/raxml-ng-adaptive"
mptp.exe_path = "./bin/mptp"
pythia.raxmlng_path = "./bin/raxml-ng"
pythia.predictor_path = "predictors/latest.pckl"
config_path = "lingdata_modeltesting_familyfull_config.json"

database.read_config(config_path)
database.compile()
df = database.data()
results_dir = os.path.join("data/results", "familyfull")
plots_dir = os.path.join("data/results", "plots")

run_raxml_ng(df)
evaluate(df)
