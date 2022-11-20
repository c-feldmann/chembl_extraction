from  typing import List, Optional
from rdkit import Chem


def gen_nonstereo_smiles(smi: str) -> Optional[str]:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def gen_inchi_key(smi: str)-> Optional[str]:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol)


def concat_unique(string_list: List[str]) -> str:
    string_set = set(string_list)
    concat_str = "_".join(sorted(string_set))
    return concat_str
