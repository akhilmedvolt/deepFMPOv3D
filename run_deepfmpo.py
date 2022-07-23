import os
from mvcommon.aws.s3 import S3
from utils import *
import glob

NCE_BUCKET = "medvolt-nce"
DATASET_NAME = "docked_smiles_qed"  # Change to cov2
MODEL_PATH = f"./data_models/chemprop/{DATASET_NAME}"


def predict_score(mol):
    convert_to_smiles(mol)
    os.system(
        f"chemprop_predict --test_path smile_df.csv --checkpoint_dir {MODEL_PATH} --preds_path ./smile_df_pred.csv"
    )
    smile_df = pd.read_csv("/smile_df_pred.csv")
    return smile_df["target"][0]


if __name__ == "__main__":
    s3 = S3()

    # Downloading & Extracting Data
    # file = f"graphinvent/{DATASET_NAME}/generated/generated_{DATASET_NAME}.csv" # Actual Location
    file = f"graphinvent/{DATASET_NAME}/generated_smiles_{DATASET_NAME}.csv"  # Temporary Location
    file_loc = f"./Data/generated_smiles_{DATASET_NAME}.csv"
    s3.download_file(
        bucket=NCE_BUCKET, key="dynamic/" + file, download_path=file_loc
    )
    csv_to_smi(file_loc)
    os.remove(file_loc)

    # Downloading & Extracting Models
    file = f"chemprop/chemprop_{DATASET_NAME}.zip"
    file_loc = f"./data_models/chemprop_{DATASET_NAME}.zip"
    s3.download_file(
        bucket=NCE_BUCKET, key="dynamic/" + file, download_path=file_loc
    )
    unzip(file_loc, f"./data_models/chemprop/{DATASET_NAME}")
    os.remove(file_loc)
    print("Files Downloaded")

    os.system("python deepFMPO.py -f ./Data/molecules.smi -l ./Data/lead.smi -o results.sdf")

    # Convert to smi files
    os.system("obabel results_uniq.sdf -osmi -O results_uniq.smi")
    os.system("obabel results.sdf -osmi -O results.smi")

    files_smi = [file for file in glob.glob("./*") if file.endswith("smi")]
    files_sdf = [file for file in glob.glob("./*") if file.endswith("sdf")]
    files = files_sdf + files_smi

    for file in files:
        s3.upload_object(
            bucket=NCE_BUCKET,
            key=f"dynamic/deepfmpo/{DATASET_NAME}/{file.split('/')[1]}",
            file_destination=f"{file.split('/')[1]}",
        )
    print("Files uploaded to S3")
