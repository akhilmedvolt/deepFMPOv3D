'''
Things to do  ?

To Calculate:
    - dock Score (-11,-9)
    - qed (0.4,0.8)
    - sascore (3,7)

Pre-Req:
    - dock Score:
        - Train the chemprop with existing data : Done
        - mol => smiles (using rdkit)
        - Now predict the score and put it in the range
    - qed:
        - direct rdkit mol input
    - sascore
        - sascore module implementation
'''

import rdkit.Chem as Chem
from rdkit.Chem import QED

from global_parameters import dock_score_upper_limit
from global_parameters import dock_score_lower_limit
from global_parameters import qed_upper_limit
from global_parameters import qed_lower_limit
from global_parameters import sas_score_upper_limit
from global_parameters import sas_score_lower_limit

from Modules.dock_predict import predict_score
from data_models.rdkit_sa.sa_score import calculate_sascore


def evaluate_chem_mol(mol):
    try:
        Chem.GetSSSR(mol)
        dock_score = predict_score(mol)
        qed = QED.qed(mol)
        sas_score = calculate_sascore(mol)
        ret_val = [
            True,
            dock_score_lower_limit < dock_score < dock_score_upper_limit,
            qed_lower_limit < qed < qed_upper_limit,
            sas_score_lower_limit < sas_score < sas_score_upper_limit
        ]
    except:
        ret_val = [False] * 4

    return ret_val
