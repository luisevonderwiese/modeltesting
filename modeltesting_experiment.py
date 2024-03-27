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
import code.distances as distances
distances.exe_path = "./bin/qdist"


def filter_data(df):
    for (i, row) in df.iterrows():
        categorical_path = row["categorical_path"]
        cd = CategoricalData.from_file(categorical_path)
        num_possible_values = [len(cd.get_possible_values(char_idx)) for char_idx in range(cd.num_chars())]
        char_ids_filtered = [char_id for (char_idx, char_id) in enumerate(cd.char_ids) if num_possible_values[char_idx] <= 2]
        if len(char_ids_filtered) < 4:
            df.drop(i, inplace=True)
            continue
        matrix_filtered = [matrix_col for (char_idx, matrix_col) in enumerate(cd.matrix) if num_possible_values[char_idx] <= 2]
        cd.char_ids = char_ids_filtered
        cd.matrix = matrix_filtered
        cd.write(categorical_path)
        cd.write_msa(row["msa_paths"]["bin"], "bin")
    return df



def run_raxml_ng(df, experiment):
    for (i, row) in df.iterrows():
        if experiment == "raxmlng_gamma":
            raxmlng.run_inference(row["msa_paths"]["bin"], "BIN+G", util.prefix(results_dir, row, experiment, "bin"))
        if experiment == "raxmlng_nogamma":
            raxmlng.run_inference(row["msa_paths"]["bin"], "BIN", util.prefix(results_dir, row, experiment, "bin"))


def run_pythia(df):
    for (i, row) in df.iterrows():
        d = "/".join(util.prefix(results_dir, row, "pythia", "bin").split("/")[:-1])
        if not os.path.isdir(d):
            os.makedirs(d)
        pythia.run_with_padding(row["msa_paths"]["bin"], util.prefix(results_dir, row, "pythia", "bin"))

def run_mptp(df):
    for (i, row) in df.iterrows():
        best_tree_path = raxmlng.best_tree_path(util.prefix(results_dir, row, "raxmlng_gamma", "bin"))
        mptp.run(best_tree_path, util.prefix(results_dir, row, "mptp_gamma", "bin"))
        best_tree_path = raxmlng.best_tree_path(util.prefix(results_dir, row, "raxmlng_nogamma", "bin"))
        mptp.run(best_tree_path, util.prefix(results_dir, row, "mptp_nogamma", "bin"))


def get_alphas(df, experiment):
    alphas = []
    for (i, row) in df.iterrows():
        alphas.append(raxmlng.alpha(util.prefix(results_dir, row, experiment, "bin")))
    return alphas

def get_final_llhs(df, experiment):
    llhs = []
    for (i, row) in df.iterrows():
        llhs.append(raxmlng.final_llh(util.prefix(results_dir, row, experiment, "bin")))
    return llhs

def get_difficulties(df):
    difficulties = []
    for (i, row) in df.iterrows():
        difficulties.append(pythia.get_difficulty(util.prefix(results_dir, row, "pythia", "bin")))
    return difficulties


def get_aics(df, experiment):
    scores = []
    for (i, row) in df.iterrows():
        scores.append(raxmlng.aic(util.prefix(results_dir, row, experiment, "bin")))
    return scores

def get_zero_base_frequencies(df, experiment):
    frequencies = []
    for (i, row) in df.iterrows():
        base_frequencies = raxmlng.base_frequencies(util.prefix(results_dir, row, experiment, "bin"))
        if len(base_frequencies) < 2:
            frequencies.append(float('nan'))
        else:
            frequencies.append(base_frequencies[0])
    return frequencies

def average_branch_length(tree_name):
    if not os.path.isfile(tree_name):
       return float("nan")
    t  =  Tree(tree_name)
    brlens = []
    for node in t.traverse("postorder"):
        brlens.append(node.dist)
    return sum(brlens)/len(brlens)

def get_num_species(df, experiment):
    return [mptp.get_num_species(util.prefix(results_dir, row, experiment, "bin")) for (_, row) in df.iterrows()]

def average_brlens(df, experiment):
    return  [average_branch_length(raxmlng.best_tree_path(util.prefix(results_dir, row, experiment, "bin"))) for (_, row) in df.iterrows()]


def calculate_relative_likelihoods(df):
    #aics_gamma <= aics_nogamma wenn +G besser
    # BIN is rl times as probable as BIN+G to minimize the information loss
    # das heißt je kleiner rl ist, um so deutlicher ist BIN+G überlegen
    names = ["AIC", "AICc", "BIC"]
    for name in names:
        relative_likelihoods = []
        for i,row in df.iterrows():
            score_gamma = row[name + "_gamma"]
            score_nogamma = row[name + "_nogamma"]
            relative_likelihoods.append(math.exp((score_gamma - score_nogamma)/2))
        df["rlh_" + name] = relative_likelihoods
    return df

def distances_glottolog(df, experiment):
    d = []
    for i, row in df.iterrows():
        glottolog_tree_path =  row["glottolog_tree_path"]
        if glottolog_tree_path == glottolog_tree_path and os.path.isfile(glottolog_tree_path):
            best_tree_path = raxmlng.best_tree_path(util.prefix(results_dir, row, experiment, "bin"))
            if os.path.isfile(best_tree_path):
                d.append(distances.gq_distance(best_tree_path, glottolog_tree_path))
            else:
                d.append(float("nan"))
        else:
            d.append(float("nan"))
    return d


def write_results_df(df):
    df["difficulty"] = get_difficulties(df)
    df["alpha"] = get_alphas(df, "raxmlng_gamma")
    df["zero_base_frequency_gamma"] = get_zero_base_frequencies(df, "raxmlng_gamma")
    df["zero_base_frequency_nogamma"] = get_zero_base_frequencies(df, "raxmlng_nogamma")
    df["final_llh_gamma"] = get_final_llhs(df, "raxmlng_gamma")
    df["final_llh_nogamma"] = get_final_llhs(df, "raxmlng_nogamma")
    scores_gamma = get_aics(df, "raxmlng_gamma")
    scores_nogamma = get_aics(df, "raxmlng_nogamma")
    df["AIC_gamma"] = [scores[0] for scores in scores_gamma]
    df["AIC_nogamma"] = [scores[0] for scores in scores_nogamma]
    df["AICc_gamma"] = [scores[1] for scores in scores_gamma]
    df["AICc_nogamma"] = [scores[1] for scores in scores_nogamma]
    df["BIC_gamma"] = [scores[2] for scores in scores_gamma]
    df["BIC_nogamma"] = [scores[2] for scores in scores_nogamma]
    avg_ml_dist_bin = []
    for i, row in df.iterrows():
        avg_ml_dist_bin.append(raxmlng.avg_ml_tree_dist(util.prefix(results_dir, row, "raxmlng_gamma", "bin")))
    df["avg_ml_tree_dist"] = avg_ml_dist_bin
    df = calculate_relative_likelihoods(df)
    df["num_species_gamma"] = get_num_species(df, "mptp_gamma")
    df["num_species_nogamma"] = get_num_species(df, "mptp_nogamma")
    df["avg_brlen_gamma"] = average_brlens(df, "raxmlng_gamma")
    df["avg_brlen_nogamma"] = average_brlens(df, "raxmlng_nogamma")
    df["gq_glottolog_gamma"] = distances_glottolog(df, "raxmlng_gamma")
    df["gq_glottolog_nogamma"] = distances_glottolog(df, "raxmlng_nogamma")
    print_df = df[["ds_id", "source", "ling_type", "family", "difficulty", "alpha", "zero_base_frequency_gamma", "zero_base_frequency_nogamma",
        "final_llh_gamma", "final_llh_nogamma", "AIC_gamma", "AIC_nogamma", "AICc_gamma", "AICc_nogamma", "BIC_gamma", "BIC_nogamma",
                    "rlh_AIC", "rlh_AICc", "rlh_BIC", "avg_ml_tree_dist", "num_species_gamma", "num_species_nogamma", "avg_brlen_gamma", "avg_brlen_nogamma",
                    "gq_glottolog_gamma", "gq_glottolog_nogamma"]]
    print_df.to_csv(os.path.join(results_dir, "raxml_pythia_results.csv"), sep = ";")

raxmlng.exe_path = "./bin/raxml-ng"
mptp.exe_path = "./bin/mptp"
pythia.raxmlng_path = "./bin/raxml-ng"
pythia.predictor_path = "predictors/latest.pckl"
config_paths = {
                "familyfull": "lingdata_modeltesting_familyfull_config.json",
                "familysplit": "lingdata_modeltesting_familysplit_config.json",
                "familyfull_filtered": "lingdata_modeltesting_familyfull_filtered_config.json",
                "familysplit_filtered": "lingdata_modeltesting_familysplit_filtered_config.json"
                }
#database.download()
for (setup, config_path) in config_paths.items():
    database.read_config(config_path)
    database.compile()
    df = database.data()
    if setup.endswith("filtered"):
        df = filter_data(df)
    results_dir = os.path.join("data/results", setup)

    #run_raxml_ng(df, "raxmlng_gamma")
    #run_raxml_ng(df, "raxmlng_nogamma")
    #run_pythia(df)
    #run_mptp(df)
    write_results_df(df)
