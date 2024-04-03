import os
import lingdata.database as database
from Bio import AlignIO
import pandas as pd

from pylotiger import get_set_partitions

def get_characters(alingment):
    characters = {}
    for record in alignment:
        record_chars = {}
        for i, char in enumerate(record.seq):
            record_chars[str(i)] = [char]
        characters[record.id] = record_chars
    return characters


config_paths = {
                "familyfull": "lingdata_modeltesting_familyfull_config.json",
                "familysplit": "lingdata_modeltesting_familysplit_config.json",
                }
#database.download()
for (setup, config_path) in config_paths.items():
    database.read_config(config_path)
    #database.compile()
    df = database.data()

    results_dir = os.path.join("data/results", setup)
    res_df = pd.DataFrame(columns = ["ds_id", "source", "ling_type", "family", \
            "entropy", "chi_distribution", "chi_percentage", "scc", \
            "entropy_vertical", "chi_distribution_vertical", "chi_percentage_vertical", "scc_vertical"])

    for i, row in df.iterrows():
        print(row["ds_id"])
        res_df.at[i, "ds_id"] = row["ds_id"]
        res_df.at[i, "source"] = row["source"]
        res_df.at[i, "ling_type"] = row["ling_type"]
        res_df.at[i, "family"] = row["family"]
        try:
            alignment = AlignIO.read(row["msa_paths"]["bin"], "phylip-relaxed")
        except:
            res_df.at[i, "entropy"] = float("nan")
            res_df.at[i, "chi_percentage"] = float("nan")
            res_df.at[i, "chi_distribution"] = float("nan")
            res_df.at[i, "scc"] = float("nan")
            res_df.at[i, "entropy_vertical"] = float("nan")
            res_df.at[i, "chi_percentage_vertical"] = float("nan")
            res_df.at[i, "chi_distribution_vertical"] = float("nan")
            res_df.at[i, "scc_vertical"] = float("nan")

            print("xxxxxxxxxxxxxxxxxxx")
            continue
        entropies = []
        chi_distributions = []
        chi_percentages = []
        sccs = []
        for record in alignment:
            with open("temp.txt", "w+") as tempfile:
                tempfile.write(str(record.seq))
            os.system("./bin/ent temp.txt > out.txt")
            with open("out.txt", "r") as outfile:
                for line in outfile.readlines():
                    if line.startswith("Entropy"):
                        try:
                            e = float(line.split(" = ")[1].split(" ")[0])
                            entropies.append(e)
                        except:
                            print("Entropy undefined")
                    if line.startswith("Chi"):
                        try:
                            c_d = float(line.split(",")[0].split(" ")[-1])
                            chi_distributions.append(c_d)
                        except:
                            print("Chi distribution undefined")
                    if line.startswith("would"):
                        try:
                            c_p = float(line.split(" percent ")[0].split(" ")[-1])
                            chi_percentages.append(c_p)
                        except:
                            print(line)
                            print("Chi percentage undefined")
                    if line.startswith("Serial correlation coefficient"):
                        try: 
                            scc = abs(float(line.split(" is ")[1].split(" ")[0]))
                            sccs.append(scc)
                        except:
                            sccs.append(1.0)
        res_df.at[i, "entropy"] = sum(entropies) / len(entropies)
        res_df.at[i, "chi_percentage"] = sum(chi_percentages) / len(chi_percentages)
        res_df.at[i, "chi_distribution"] = sum(chi_distributions) / len(chi_distributions)
        res_df.at[i, "scc"] = sum(sccs) / len(sccs)

        entropies = []
        chi_distributions = []
        chi_percentages = []
        sccs = []

        for k  in range(alignment.get_alignment_length()):
            site = alignment[:, k]
            with open("temp.txt", "w+") as tempfile:
                tempfile.write(str(site))
            os.system("./bin/ent temp.txt > out.txt")
            with open("out.txt", "r") as outfile:
                for line in outfile.readlines():
                    if line.startswith("Entropy"):
                        try:
                            e = float(line.split(" = ")[1].split(" ")[0])
                            entropies.append(e)
                        except:
                            print("Entropy undefined")
                    if line.startswith("Chi"):
                        try:
                            c_d = float(line.split(",")[0].split(" ")[-1])
                            chi_distributions.append(c_d)
                        except:
                            print("Chi distribution undefined")
                    if line.startswith("would"):
                        try:
                            c_p = float(line.split(" percent ")[0].split(" ")[-1])
                            chi_percentages.append(c_p)
                        except:
                            print("Chi percentage undefined")
                    if line.startswith("Serial correlation coefficient"):
                        try:
                            scc = abs(float(line.split(" is ")[1].split(" ")[0]))
                            sccs.append(scc)
                        except:
                            sccs.append(1.0)
        res_df.at[i, "entropy_vertical"] = sum(entropies) / len(entropies)
        res_df.at[i, "chi_percentage_vertical"] = sum(chi_percentages) / len(chi_percentages)
        res_df.at[i, "chi_distribution_vertical"] = sum(chi_distributions) / len(chi_distributions)
        res_df.at[i, "scc_vertical"] = sum(sccs) / len(sccs)


    res_df.to_csv(os.path.join(results_dir, "ent.csv"), sep = ";")
