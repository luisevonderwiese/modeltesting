import os
import lingdata.database as database
from Bio import AlignIO
from nexus import NexusReader
from phylogemetric import DeltaScoreMetric
from phylogemetric import QResidualMetric


config_paths = {
                "familyfull": "lingdata_modeltesting_familyfull_config.json",
                "familysplit": "lingdata_modeltesting_familysplit_config.json",
                }
#database.download()
for (setup, config_path) in config_paths.items():
    database.read_config(config_path)
    #database.compile()
    df = database.data()
    if setup.endswith("filtered"):
        df = filter_data(df)
    results_dir = os.path.join("data/results", setup)
    res_df = pd.DataFrame(columns = ["ds_id", "source", "ling_type", "family", "q", "delta"])
    deltas = []
    qs = []
    for i, row in df.iterrows():
        res_df.at[i, "ds_id"] = row["ds_id"]
        res_df.at[i, "source"] = row["source"]
        res_df.at[i, "ling_type"] = row["ling_type"]
        res_df.at[i, "family"] = row["family"]
        alignment = AlignIO.read(row["msa_paths"]["bin"], "phylip-relaxed")
        with open("temp.nex", "w+") as outfile:
            AlignIO.write(alignment, outfile, "nexus")
            matrix = NexusReader("temp.nex").data.matrix
            res_df.at[i, "q"] = QResidualMetric(matrix)
            res_df.at[i, "delta"] = DeltaScoreMetric(matrix)
    res_df.to_csv(os.path.join(results_dir, "simon_metrics.csv"), sep = ";")
