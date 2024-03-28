import os
import lingdata.database as database
from Bio import AlignIO
import pandas as pd

from pylotiger import get_set_partitions, get_rates, corrected_pas

def get_characters(alingment):
    characters = {}
    taxa_names = [record.id for record in alignment]
    for c in range(alignment.get_alignment_length()):
        current = {}
        for r, taxa_name in enumerate(taxa_names):
            current[taxa_name] = alignment[r,c]
        characters[str(c)] = current
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
    res_df = pd.DataFrame(columns = ["ds_id", "source", "ling_type", "family", "tiger", "tiger_corrected"])

    for i, row in df.iterrows():
        res_df.at[i, "ds_id"] = row["ds_id"]
        res_df.at[i, "source"] = row["source"]
        res_df.at[i, "ling_type"] = row["ling_type"]
        res_df.at[i, "family"] = row["family"]
        try:
            alignment = AlignIO.read(row["msa_paths"]["bin"], "phylip-relaxed")
        except:
            res_df.at[i, "tiger"] = float("nan")
            continue
        characters = get_characters(alignment)
        taxa = [record.id for record in alignment]
        set_partitions = get_set_partitions(characters, taxa)
        rates = get_rates(set_partitions)
        res_df.at[i, "tiger"] = sum(rates.values()) / len(rates)
        corrected_rates = []
        for i1, p1 in set_partitions.items():
            for i2, p2 in set_partitions.items():
                if i1 == i2:
                    continue
                pas = corrected_pas(p1, p2)
                if pas is None:
                    continue
                else:
                    corrected_rates.append(pas)
        if len(corrected_rates) == 0:
            res_df.at[i, "tiger_corrected"] = float("nan")
        else:
            res_df.at[i, "tiger_corrected"] = sum(corrected_rates) / len(corrected_rates)
    res_df.to_csv(os.path.join(results_dir, "tiger.csv"), sep = ";")
