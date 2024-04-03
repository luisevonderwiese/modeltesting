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
    res_df = pd.DataFrame(columns = ["ds_id", "source", "ling_type", "family", "entropy", "chi_percentage", "scc"])

    for i, row in df.iterrows():
        res_df.at[i, "ds_id"] = row["ds_id"]
        res_df.at[i, "source"] = row["source"]
        res_df.at[i, "ling_type"] = row["ling_type"]
        res_df.at[i, "family"] = row["family"]
        try:
            alignment = AlignIO.read(row["msa_paths"]["bin"], "phylip-relaxed")
        except:
            res_df.at[i, "entropy"] = float("nan")
            res_df.at[i, "chi_percentage"] = float("nan")
            res_df.at[i, "scc"] = float("nan")
            continue
        for record in alignment:
            with open("temp.txt", "w+") as tempfile:
                tempfile.write(str(record.seq))
            #print(open("temp.txt", "r").readlines())
            os.system("./ent temp.txt > out.txt")
            with open("out.txt", "r") as outfile:
                for line in outfile.readlines():
                    if line.startswith("Entropy"):
                        entropies.append(float(line.split(" = ")[1].split(" ")[0]))
                    if line.startswith("would"):
                        chi_percentages.append(float(line.split(" percent ")[0].split(" ")[-1]))
                    if line.startswith("Serial correlation coefficient"):
                        sccs.append(float(line.split(" is ")[1].split(" ")[0]))
        res_df.at[i, "entropy"] = sum(entropies) / len(entropies)
        res_df.at[i, "chi_percentage"] = sum(chi_percentages) / len(chi_percentages)
        res_df.at[i, "scc"] = sum(sccs) / len(scc)

    res_df.to_csv(os.path.join(results_dir, "ent.csv"), sep = ";")
