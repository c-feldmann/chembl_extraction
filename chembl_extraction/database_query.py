import pandas as pd
import sqlite3


POTENCY_DATA_QUERY = """
SELECT DISTINCT 
 compound_structures.canonical_smiles,
 molecule_dictionary.chembl_id as chembl_cid,
 compound_records.compound_key,
 component_sequences.accession as uniprot_id,
 component_sequences.db_source, 
 target_dictionary.chembl_id as chembl_tid,  
 target_dictionary.pref_name, 
 target_dictionary.organism,
 variant_sequences.mutation,
 variant_sequences.isoform,
 activities.activity_id,
 activities.standard_value,
 activities.standard_units,
 activities.standard_relation,
 activities.standard_type,
 activities.activity_comment,
 activities.data_validity_comment,
 assays.chembl_id as chembl_aid,
 bioassay_ontology.label as bao_format,
 docs.chembl_id as chembl_docid,
 docs.doc_type
FROM 
 activities, 
 assays,
 bioassay_ontology,
 compound_records,
 compound_structures,
 component_sequences,
 docs,
 molecule_dictionary,
 molecule_hierarchy,
 target_dictionary, 
 target_components
LEFT JOIN variant_sequences
ON assays.variant_id = variant_sequences.variant_id
WHERE 
 activities.assay_id = assays.assay_id AND
 activities.doc_id = docs.doc_id AND
 activities.molregno = molecule_dictionary.molregno AND
 activities.molregno = molecule_hierarchy.molregno AND
 activities.record_id = compound_records.record_id AND
 assays.tid = target_dictionary.tid AND
 assays.bao_format = bioassay_ontology.bao_id AND
 molecule_hierarchy.parent_molregno = compound_structures.molregno AND
 target_dictionary.tid = target_components.tid AND
 target_components.component_id = component_sequences.component_id AND
 assays.relationship_type = 'D' AND 
 assays.confidence_score = 9 AND 
 target_dictionary.target_type = 'SINGLE PROTEIN' AND
 (activities.standard_type = 'Ki' OR 
  activities.standard_type = 'IC50' OR
  activities.standard_type = 'Kd')
"""


def extract_potency_data(connection: sqlite3.Connection) -> pd.DataFrame:
    chembl_df: pd.DataFrame = pd.read_sql(POTENCY_DATA_QUERY, connection)
    # Curating Drug matrix data
    a_comment = 'Not Active (inhibition < 50% @ 10 uM and thus dose-reponse curve not measured)'
    chembl_df.loc[chembl_df.activity_comment == a_comment, 'standard_value'] = 10000  # 10µM = 10,000 nM
    chembl_df.loc[chembl_df.activity_comment == a_comment, 'standard_relation'] = '>'
    chembl_df.loc[chembl_df.activity_comment == a_comment, 'standard_units'] = 'nM'
    chembl_df = chembl_df.query("standard_units == 'nM'")
    chembl_df = chembl_df.query("~standard_relation.isna()")
    chembl_df = chembl_df.query("~standard_value.isna()")
    chembl_df = chembl_df.query("standard_relation != '~'")
    chembl_df.loc[chembl_df.standard_relation == ">=", "standard_relation"] = ">"
    chembl_df.loc[chembl_df.standard_relation == "<=", "standard_relation"] = "<"
    chembl_df.loc[chembl_df.standard_relation == ">>", "standard_relation"] = ">"
    chembl_df.loc[chembl_df.standard_relation == "<<", "standard_relation"] = "<"

    return chembl_df


if __name__ == "__main__":
    con = sqlite3.connect("../data/chem")
