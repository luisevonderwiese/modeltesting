import os
import lingdata.database as database
from Bio import AlignIO
from phylogemetric import DeltaScoreMetric
from phylogemetric import QResidualMetric
import pandas as pd

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
    res_df = pd.DataFrame(columns = ["ds_id", "source", "ling_type", "family", "q", "delta"])

    for i, row in df.iterrows():
        res_df.at[i, "ds_id"] = row["ds_id"]
        res_df.at[i, "source"] = row["source"]
        res_df.at[i, "ling_type"] = row["ling_type"]
        res_df.at[i, "family"] = row["family"]
        try:
            alignment = AlignIO.read(row["msa_paths"]["bin"], "phylip-relaxed")
        except:
            res_df.at[i, "q"] = float("nan")
            res_df.at[i, "delta"] = float("nan")
            continue

        matrix =  {}
        for record in alignment:
            matrix[record.id] = record.seq
        try:
            q_values = QResidualMetric(matrix).score(workers=4)
            res_df.at[i, "q"] = sum(q_values.values()) / len(q_values)
        except:
            res_df.at[i, "q"] = float('nan')
        try:
            delta_scores = DeltaScoreMetric(matrix).score(workers=4)
            res_df.at[i, "delta"] = sum(delta_scores.values()) /  len(delta_scores)
        except:
            res_df.at[i, "delta"] = float('nan')
    res_df.to_csv(os.path.join(results_dir, "simon_metrics.csv"), sep = ";")
