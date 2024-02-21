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



raxmlng.exe_path = "./bin/raxml-ng"
confing_path = "lingdata_modeltesting_familyfull_partitioning.json"
database.read_config(config_path)
database.compile()
df = database.data()
results_dir = os.path.join("data/results", setup)

run_raxml_ng(df, "raxmlng_gamma")
