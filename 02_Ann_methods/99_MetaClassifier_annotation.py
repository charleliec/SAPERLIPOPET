import os
import numpy as np
import pandas as pd
import re


dfs = []
for path, model in zip(snakemake.input["metadata_csv"], snakemake.params["models"]):
    df = pd.read_csv(path, index_col=0)

    # rename only the varying columns
    df = df.rename(
        columns={
            "predictions": f"predictions_{model}",
            "mapqc_score": f"mapqc_score_{model}",
            "mapqc_score_binary": f"binary_mapqc_score_{model}",
        }
    )
    dfs.append(df)

# concat on index
merged = pd.concat(dfs, axis=1)

# if common columns are duplicated, keep only one copy
merged = merged.loc[:, ~merged.columns.duplicated()]


pred_cols = [c for c in merged.columns if c.startswith("predictions_")]
def row_consensus(row):
    counts = row[pred_cols].value_counts()
    total = counts.sum()
    top1 = counts.iloc[0]
    top2 = counts.iloc[1] if len(counts) > 1 else 0
    return pd.Series({
        "final_prediction": counts.index[0],
        "pred_label_proportion": top1 / total,
        "proportion_difference": (top1 - top2) / total
    })
merged = merged.join(merged.apply(row_consensus, axis=1))


score_cols = [c for c in df.columns if c.startswith("mapqc_score_")]
def row_scores(row):
    vals = row[score_cols].astype(float)
    return pd.Series({
        "mapqc_median": vals.median(),
        "mapqc_q3": vals.quantile(0.75),
        "mapqc_max": vals.max()
    })
merged = merged.join(merged.apply(row_scores, axis=1))

map_qc_cols = ["mapqc_median", "mapqc_q3", "mapqc_max"]

for col in map_qc_cols:
    new_c = "binary_" + col
    merged[new_c] = np.where(merged[col].isna(), "NaN",
                  np.where(merged[col] > 2, ">2", "<=2"))



merged.to_csv(snakemake.output["meta_metadata"])
print("MetaClassifier prediction saved in : ", snakemake.output["meta_metadata"])