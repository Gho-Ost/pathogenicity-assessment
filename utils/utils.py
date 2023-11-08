"""
Module for custom utility functions.
"""

import pandas as pd
from pathlib import Path

def get_data(sample_folder_path, file_types=None):
    """
    Usage example:
    
    from utils.utils import get_data
    
    EE_015 = get_data("data/EE_015/", ["default", "csq"])
    EE_050 = get_data("data/EE_050/", ["default", "csq"])
    EE_069 = get_data("data/EE_069/", ["default", "csq"])
    df = pd.concat([EE_015, EE_050, EE_069], ignore_index=True, axis=0)
    """
    types = ["default", "csq", "genotype"]
    samples = ["EE_015", "EE_050", "EE_069"]

    if file_types is None:
        raise Exception(f"No file types provided") 

    for f_type in file_types:
        if f_type not in types:
            raise Exception(f"File type {f_type} incorrect") 
        
    folder_path = Path(sample_folder_path)
    sample = sample_folder_path.rstrip("/\\")[-6:]
    
    if sample not in samples:
        raise Exception(f"Sample {sample} incorrect")

    if len(file_types) == 1:
        return pd.read_csv(folder_path / f"{sample}_{file_types[0]}.csv.gz", sep=";", compression="gzip").drop("Unnamed: 0", axis=1)
    
    df = pd.read_csv(folder_path / f"{sample}_{file_types[0]}.csv.gz", sep=";", compression="gzip").drop("Unnamed: 0", axis=1)
    for f_type in file_types[1:]:
        df2 = pd.read_csv(folder_path / f"{sample}_{f_type}.csv.gz", sep=";", compression="gzip").drop("Unnamed: 0", axis=1)
        df = pd.concat([df, df2], axis=1, ignore_index=False)
        
    return df

        
    
