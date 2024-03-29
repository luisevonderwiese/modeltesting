import os
from ete3 import Tree
import code.distances as distances

exe_path = ""


def best_tree_path(prefix):
    return prefix + ".raxml.bestTree"

def ml_trees_path(prefix):
    return prefix + ".raxml.mlTrees"

def siterate_lhs_path(prefix):
    return prefix + ".raxml.siterateLH"

def siterate_lhs(prefix):
    path = siterate_lhs_path(prefix)
    if not os.path.isfile(path):
        return []
    with open(path, "r") as f:
        lines = f.readlines()
    siteratelhs = []
    for line in lines:
        siteratelhs.append([float(el) for el in line.split(" ")[:-1]])
    return siteratelhs

def avg_ml_tree_dist(prefix):
    path = ml_trees_path(prefix)
    if not os.path.isfile(path):
        return float('nan')
    trees = [Tree(l) for l in open(path, "r").readlines()]
    dists = []
    for (i, tree1) in enumerate(trees):
        for j in range(i+1, len(trees)):
            tree2 = trees[j]
            dists.append(distances.rf_distance(tree1, tree2))
    return sum(dists) / len(dists)

def bootstrap_trees_path(prefix):
    return prefix + ".raxml.bootstraps"

def support_tree_path(prefix):
    return prefix + ".raxml.support"

def support_values(prefix):
    if not os.path.isfile(support_tree_path(prefix)):
        return []
    support_tree = Tree(support_tree_path(prefix))
    supports = [node.dist / 100 for node in support_tree.traverse("postorder")]
    return supports


def consensus_tree_path(prefix):
    return prefix + ".raxml.consensusTreeMR"

def alpha(prefix):
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("   Rate heterogeneity:"):
            return float(line.split(",  ")[1].split(" ")[1])
    return float('nan')

def weights_rates(prefix):
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("   Rate heterogeneity:"):
            string = line.split(",  ")[2].split(": ")[1]
            parts = string.split(" ")[:4]
            res = []
            for part in parts:
                part = part.strip("()")
                res.append([float(inner_part) for inner_part in part.split(",")])
            return res
    return []

def final_llh(prefix):
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Final LogLikelihood: "):
            return float(line.split(": ")[1])
    return float('nan')

def aic(prefix):
    logpath = prefix + ".raxml.log"
    if not os.path.isfile(logpath):
        return [float('nan'), float('nan'), float('nan')]
    with open(logpath, "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("AIC"):
            parts = line.split(" / ")
            scores = []
            for part in parts:
                scores.append(float(part.split(" ")[2]))
            return scores
    return [float('nan'), float('nan'), float('nan')]


def get_partitioning_alphas(prefix):
    file_name = prefix + ".raxml.bestModel"
    if not os.path.isfile(file_name):
        return {}
    with open(file_name, "r") as logfile:
        lines = logfile.readlines()
    alphas = {}
    for line in lines:
        idx = int(line.split(" ")[1].split("p")[1])
        alpha = float(line.split("{")[2].split("}")[0])
        alphas[idx] = alpha
    return alphas



def base_frequencies(prefix):
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("   Base frequencies"):
            parts = line.split(": ")[1].split(" ")[:-1]
            return [float(part) for part in parts]
    return []

def consense_tree(prefixes, prefix, args = ""):
    if not os.path.isfile(consensus_tree_path(prefix)):
        args = args + " --redo"
    with open("trees.nw", "w+") as outfile:
        for p in prefixes:
            tree_path = best_tree_path(p)
            if os.path.isfile(tree_path):
                with open(tree_path, "r") as infile:
                    outfile.write(infile.read())
    command = exe_path
    command += " --consense"
    command += " --tree trees.nw "
    command += " --prefix " + prefix
    os.system(command)
    os.remove("trees.nw")


def run_inference(msa_path, model, prefix, args = ""):
    if exe_path == "":
        print("Please specify raxmlng.exe_path")
        return
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(best_tree_path(prefix)):
        args = args + " --redo"
    command = exe_path
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2 --force model_lh_impr"
    command += " " + args
    os.system(command)

def run_bootstrap(msa_path, model, prefix, args = ""):
    if exe_path == "":
        print("Please specify raxmlng.exe_path")
        return
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(bootstrap_trees_path(prefix)):
        args = args + " --redo"
    command = exe_path
    command += " --bootstrap"
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2"
    command += " " + args
    os.system(command)


def run_support(inference_prefix, bootstrap_prefix, prefix, args = ""):
    if exe_path == "":
        print("Please specify raxmlng.exe_path")
        return
    if os.path.isfile(prefix + ".log"):
        print("Run files for " + prefix + " present")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(support_tree_path(prefix)):
        args = args + " --redo"
    command = exe_path
    command += " --support"
    command += " --tree " + best_tree_path(inference_prefix)
    command += " --bs-trees " + bootstrap_trees_path(bootstrap_prefix)
    command += " --prefix " + prefix
    command += " --threads 2 "
    command += " " + args
    os.system(command)


def run_bootstrap(msa_path, model, prefix, args = ""):
    if exe_path == "":
        print("Please specify raxmlng.exe_path")
        return
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(bootstrap_trees_path(prefix)):
        args = args + " --redo"
    command = exe_path
    command += " --bootstrap"
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2"
    command += " " + args
    os.system(command)


def run_inference_adaptive(msa_path, model, prefix, args = ""):
    if adaptive_exe_path == "":
        print("Please specify raxmlng.exe_path")
        return
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(best_tree_path(prefix)):
        args = args + " --redo"
    command = adaptive_exe_path
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --search --thread 1 --simd none --site-repeats off --pat-comp off"
    command += " " + args
    os.system(command)
