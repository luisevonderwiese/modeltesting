import os

def prefix(results_dir, row, experiment, run):
    return os.path.join(results_dir, experiment, "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]]), run)
