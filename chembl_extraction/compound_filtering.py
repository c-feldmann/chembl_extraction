import io
import pandas as pd
import os
import subprocess
from chembl_extraction.utils import concat_unique, gen_inchi_key, gen_nonstereo_smiles

from compchemkit.filtering import PainsFilter
from chembl_extraction.variables import LILLY_MEDCHEM_PATH
from chembl_extraction.variables import LILLY_MEDCHEM_WORKING_PATH
from chembl_extraction.variables import AGGREGATOR_PATH


def gen_compound_df(chembl_df: pd.DataFrame) -> pd.DataFrame:
    cpd_df = chembl_df.groupby("canonical_smiles")[["chembl_cid"]].agg(concat_unique)
    cpd_df.reset_index(inplace=True)
    cpd_df["nonstereo_smiles"] = cpd_df.canonical_smiles.apply(gen_nonstereo_smiles)
    cpd_df["inchi_key"] = cpd_df.canonical_smiles.apply(gen_inchi_key)
    cpd_df["inchi_key_head"] = cpd_df["inchi_key"].apply(lambda x: x.split("-")[0])

    pains_filter = PainsFilter()
    cpd_df["is_PAINS"] = pains_filter.check_smiles_list(cpd_df.canonical_smiles.tolist())

    aggregator_df = pd.read_csv(AGGREGATOR_PATH, sep=" ", names=["smiles", "zinc_id"])
    aggregator_df["inchi_key"] = aggregator_df.smiles.apply(gen_inchi_key)
    aggregator_df["inchi_key_head"] = aggregator_df["inchi_key"].apply(lambda x: x.split("-")[0])
    cpd_df["is_aggregator"] = cpd_df["inchi_key_head"].isin(aggregator_df["inchi_key_head"])

    cpd_df.to_csv(
        os.path.join(LILLY_MEDCHEM_WORKING_PATH, "unique_cpds.smi"),
        sep=' ',
        index=False,
        columns=["nonstereo_smiles", "chembl_cid"])

    curr_dir = os.path.abspath(".")
    try:
        os.chdir(LILLY_MEDCHEM_WORKING_PATH)
        out = subprocess.run(["ruby", LILLY_MEDCHEM_PATH, "unique_cpds.smi"], capture_output=True)
    finally:
        os.chdir(curr_dir)
    lilly_output = io.StringIO(out.stdout.decode())
    lilly_df = pd.read_csv(lilly_output, sep=" ", names=["lilly_smiles", "chembl_cid", "I", "dont", " know"])

    cpd_df["is_lilly_flagged"] = True
    cpd_df.loc[lambda row: row["chembl_cid"].isin(lilly_df.chembl_cid), "is_lilly_flagged"] = False

    cpd_df["is_interference"] = False
    cpd_df.loc[cpd_df.is_PAINS, "is_interference"] = True
    cpd_df.loc[cpd_df.is_lilly_flagged, "is_interference"] = True
    cpd_df.loc[cpd_df.is_aggregator, "is_interference"] = True
    cpd_df.drop(columns=["inchi_key_head"], inplace=True)
    return cpd_df
