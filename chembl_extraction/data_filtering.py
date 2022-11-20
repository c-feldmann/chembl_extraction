import pandas as pd


def filter_med_conf(chembl_df: pd.DataFrame) -> pd.DataFrame:
    chembl_df = chembl_df.query("data_validity_comment != 'Potential transcription error'")
    chembl_df = chembl_df.query("data_validity_comment != 'Potential author error'")
    return chembl_df

def filter_high_conf(chembl_df: pd.DataFrame) -> pd.DataFrame:
    chembl_df = chembl_df.query("doc_type != 'PATENT'")
    chembl_df = chembl_df.query("chembl_docid != 'CHEMBL1201862'")  # Remove Pubchem
    chembl_df = chembl_df.query("mutation.isna()")
    chembl_df = chembl_df.query("bao_format == 'single protein format'")
    return chembl_df
