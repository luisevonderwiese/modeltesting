import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import re
from collections import Counter
from tabulate import tabulate

from lingdata.categorical import CategoricalData
from lingdata import database



def get_swadesh_ratio(categorical_path, swadesh100, swadesh207):
    raw_char_ids = CategoricalData.from_file(categorical_path).char_ids
    char_ids = set([re.sub('[^a-zA-Z]', '', c).lower() for c in raw_char_ids])
    if len(char_ids) != len(raw_char_ids):
        return float('nan')
    num_swadesh_words = len([char_id for char_id in char_ids if char_id in swadesh100 or char_id in swadesh207])
    return num_swadesh_words / len(char_ids)



def get_bins(arr, nbins):
    min_val = min(arr)
    max_val = max(arr)
    step = (max_val - min_val) / (nbins - 1)
    return [min_val + i * step for i in range(nbins)]



def AIC_analysis(df):
    gamma = 0
    no_gamma = 0
    gamma_rlhs = []
    nogamma_rlhs = []
    for i, row in df.iterrows():
        if row["AIC_gamma"] != row["AIC_gamma"] or row["AIC_nogamma"] != row["AIC_nogamma"]:
            continue
        if row["AIC_gamma"] <= row["AIC_nogamma"]:
            gamma_rlhs.append(math.exp((row["AIC_gamma"]  - row["AIC_nogamma"])/2))
            gamma += 1
        else:
            nogamma_rlhs.append(math.exp((row["AIC_nogamma"]  - row["AIC_gamma"])/2))
            no_gamma += 1
    print("Number of datasets for which BIN has the lower AIC score: " + str(no_gamma))
    print("Average relative likelihood: " + str(sum(nogamma_rlhs) / len(nogamma_rlhs)))
    print("Number of datasets for which BIN+G has the lower AIC score: " + str(gamma))
    print("Average relative likelihood: " + str(sum(gamma_rlhs) / len(gamma_rlhs)))


    # names = ["AIC"]
    # for name in names:
    #     fig, ax = plt.subplots()
    #     ax.scatter(df["alpha"], df["rlh_" + name], s=2)
    #     ax.set_yscale('log')
    #     plt.xlabel('alpha')
    #     plt.ylabel('relative likelihood of BIN')
    #     plt.savefig(os.path.join(plots_dir, "relative_lhs_" + name + '.png'))
    #     plt.clf()


# def alpha_analysis(df):
#     plt.hist(df["alpha"], bins = get_bins(df["alpha"], 50))
#     plt.xlabel('alpha')
#     plt.ylabel('Number of datasets')
#     plt.savefig(os.path.join(plots_dir, 'alphas.png'))
#     plt.clf()


def familysplit_analysis(df):
    all_families = set()
    high_alpha_families = {}
    low_alpha_families = {}
    for i, row in df.iterrows():
        sub_families = row["sub_families"]
        if (len(sub_families) != 1):
            continue
        all_families.update(sub_families)
        for f in sub_families:
            if row["alpha"] < 40:
                if f in low_alpha_families:
                    low_alpha_families[f] += 1
                else:
                    low_alpha_families[f] = 1
            else:
                if f in high_alpha_families:
                    high_alpha_families[f] += 1
                else:
                    high_alpha_families[f] = 1

    columns = ["family", "num_high_alpha", "num_low_alpha"]
    matrix = []
    for f in all_families:
        row = [f]
        if f in high_alpha_families:
            row.append(high_alpha_families[f])
        else:
            row.append(0)
        if f in low_alpha_families:
            row.append(low_alpha_families[f])
        else:
            row.append(0)
        matrix.append(row)
    new_df = pd.DataFrame(matrix, columns = columns)
    num_mixed = len(set(high_alpha_families.keys()).intersection(set(low_alpha_families.keys())))
    print(new_df)

    low_cnt = 0
    high_cnt = 0
    both_more_high = 0
    both_more_low = 0
    both_equal = 0
    for i, row in new_df.iterrows():
        if row["num_high_alpha"] != 0:
            if row["num_low_alpha"] != 0:
                if row["num_high_alpha"] > row["num_low_alpha"]:
                    both_more_high += 1
                elif row["num_high_alpha"] < row["num_low_alpha"]:
                    both_more_low += 1
                else:
                    both_equal += 1
            else:
                high_cnt += 1
        else:
            low_cnt += 1
    print("High alpha only: " +  str(high_cnt))
    print("Low alpha only: " +  str(low_cnt))
    print("Both, more high: " +  str(both_more_high))
    print("Both, more low: " +  str(both_more_low))
    print("Both, equal: " +  str(both_equal))

def familyfull_analysis(df):
    num_families_low = []
    num_families_high = []
    for i, row in df.iterrows():
        sub_families = row["sub_families"]
        if row["alpha"] < 40:
            num_families_low.append(len(sub_families))
        else:
            num_families_high.append(len(sub_families))
    print("Number of datasets with low alpha (high rate heterogenity) containing x subfamilies")
    print(Counter(num_families_low))
    print("Number of datasets with high alpha (low rate heterogenity) containing x subfamilies")
    print(Counter(num_families_high))


def alpha_correlation(columns, familyfull_df, familysplit_df):
    r = [[column, familyfull_df[column].corr(familyfull_df["alpha"]), familysplit_df[column].corr(familysplit_df["alpha"])] for column in columns]
    print(tabulate(r, tablefmt="pipe", floatfmt=".3f", headers = ["column", "alpha corr full", "alpha corr split"]))



def add_results(df):
    results_df = pd.read_csv(os.path.join(results_dir, "raxml_pythia_results.csv"), sep = ";")
    df = pd.merge(df, results_df, how = 'left', left_on=["ds_id", "source", "ling_type", "family"], right_on = ["ds_id", "source", "ling_type", "family"])
    #ebg_df = pd.read_csv(os.path.join(results_dir, "EBG_features.csv"), sep = ";")
    #df = pd.merge(df, ebg_df, how = 'left', left_on=["ds_id", "source", "ling_type", "family"], right_on = ["ds_id", "source", "ling_type", "family"])
    return df



with open("word_lists/swadesh100.txt") as s100_file:
   swadesh100 = s100_file.read().split("\n")[:-1]

with open("word_lists/swadesh207.txt") as s207_file:
   swadesh207 = s207_file.read().split("\n")[:-1]



config_paths = {"familyfull": "lingdata_modeltesting_familyfull_config.json", "familysplit": "lingdata_modeltesting_familyfull_config.json"}

dfs = {}
for (setup, config_path) in config_paths.items():
    database.read_config(config_path)
    df = database.data()
    results_dir = os.path.join("data/results", setup)
    plots_dir = os.path.join(results_dir, "plots")

    df = add_results(df)
    df["swadesh_ratio"] = [get_swadesh_ratio(row["categorical_path"], swadesh100, swadesh207) for i, row in df.iterrows()]

    df["num_species_ratio_gamma"] = df["num_species_gamma"] / df["num_taxa"]
    df["num_species_ratio_nogamma"] = df["num_species_nogamma"] / df["num_taxa"]

    dfs[setup] = df

familyfull_analysis(dfs["familyfull"])
familysplit_analysis(dfs["familysplit"])

columns = [
                "num_taxa",
                "num_chars",
                "multistate_ratio",
                "max_values",
                "difficulty",
                "AIC_gamma",
                "swadesh_ratio",
                "sites_per_char",
               # "mean_substitution_frequency",
               # "mean_norm_rf_distance",
               # "mean_parsimony_support",
               # "mean_parsimony_bootstrap_support",
                "avg_brlen_gamma",
                "avg_brlen_nogamma",
                "num_species_gamma",
                "num_species_nogamma",
                "num_species_ratio_gamma",
                "num_species_ratio_nogamma",
                ]
alpha_correlation(columns, dfs["familyfull"], dfs["familysplit"])
print("full:")
AIC_analysis(dfs["familyfull"])
print("split:")
AIC_analysis(dfs["familysplit"])
