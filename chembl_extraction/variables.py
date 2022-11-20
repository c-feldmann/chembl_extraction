import os

LILLY_MEDCHEM_PATH = "/home/bob/software/Lilly-Medchem-Rules/Lilly_Medchem_Rules.rb"
assert os.path.isfile(LILLY_MEDCHEM_PATH)

# Can't touch this! (Except you know what you are doing.)
LILLY_MEDCHEM_WORKING_PATH = os.path.abspath(os.path.join(__file__, "../../data/tmp"))
assert os.path.isdir(LILLY_MEDCHEM_WORKING_PATH), LILLY_MEDCHEM_WORKING_PATH

# Can't touch this! (Except you know what you are doing.)
AGGREGATOR_PATH = os.path.abspath(os.path.join(__file__, "../../data/aggregators.smi"))
assert os.path.isfile(AGGREGATOR_PATH), AGGREGATOR_PATH
