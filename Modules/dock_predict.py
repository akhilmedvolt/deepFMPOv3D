import pandas as pd
from rdkit import Chem
import os

MODEL_PATH = "./data_models/chemprop/delaney_checkpoints/fold_0/model_0/model.pt"


def convert_to_smiles(mol):
    smile = Chem.MolToSmiles(mol)
    smile_df = pd.DataFrame({"smiles": smile})
    smile_df.to_csv("./smile_df.csv")


def predict_score(mol):
    convert_to_smiles(mol)
    os.system(f"chemprop_predict --test_path smile_df.csv --checkpoint_path {MODEL_PATH} --preds_path "
              f"./smile_df_pred.csv")
    smile_df = pd.read_csv("/smile_df_pred.csv")
    return smile_df["logSolubility"][0]
