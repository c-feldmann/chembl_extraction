import pandas as pd
from typing import List, Optional


def aggregate_data(activity_df: pd.DataFrame, merge_columns: Optional[List[str]]=None) -> pd.DataFrame:
    if merge_columns is None:
        merge_columns = ["nonstereo_smiles", "uniprot_id", "standard_type", "standard_relation"]
    activity_agg_df = activity_df.groupby(["nonstereo_smiles", "uniprot_id", "standard_type", "standard_relation"]).agg(
        {"pPot": ["min", "mean", "max"]})
    activity_agg_df = activity_agg_df.reset_index("standard_relation")

    higher_pot_than = activity_agg_df.loc[lambda row: row[("standard_relation", "")] == "<"]
    higher_pot_than = higher_pot_than.drop(columns=[("pPot", "min"), ("pPot", "mean")])
    higher_pot_than.columns = [col[0] for col in higher_pot_than.columns]
    higher_pot_than.rename(columns={"pPot": "higher_than"}, inplace=True)
    higher_pot_than = higher_pot_than.loc[lambda row: row["higher_than"] > 3.5]

    lower_pot_than = activity_agg_df.loc[lambda row: row[("standard_relation", "")] == ">"]
    lower_pot_than = lower_pot_than.drop(columns=[("pPot", "max"), ("pPot", "mean")])
    lower_pot_than.columns = [col[0] for col in lower_pot_than.columns]
    lower_pot_than.rename(columns={"pPot": "lower_than"}, inplace=True)
    lower_pot_than = lower_pot_than.loc[lambda row: row["lower_than"] > 3.5]

    equal_pot_estimate = activity_agg_df.loc[lambda row: row[("standard_relation", "")] == "="]
    equal_pot_estimate = equal_pot_estimate.loc[lambda row: row[("pPot", "max")] < (row[("pPot", "min")] + 1)]
    equal_pot_estimate = equal_pot_estimate.drop(columns=[("pPot", "max"), ("pPot", "min")])
    equal_pot_estimate.columns = [col[0] for col in equal_pot_estimate.columns]
    equal_pot_estimate = equal_pot_estimate.loc[lambda row: row["pPot"] > 3.5]
    combined_df = lower_pot_than[["lower_than"]].merge(
        equal_pot_estimate["pPot"], left_index=True, right_index=True, how="outer")
    combined_df = combined_df.merge(
        higher_pot_than["higher_than"], left_index=True, right_index=True, how="outer")
    return combined_df
