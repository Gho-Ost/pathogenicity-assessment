"""
Module for custom utility functions.
"""

import pandas as pd
from pathlib import Path


def get_dataset(data_folder, samples, file_types, option_csq, options_genotype):
    """
    Returns cleaned up dataset from chosen samples and file types.

    Parameters
    ----------
    data_folder : path to data folder e.g. "./data"

    samples : samples to be included in the dataset e.g. ["EE_015", "EE_050"] possible options: ["EE_015", "EE_050", "EE_069"]

    file_types : file types to be included in the dataset e.g. ["csq"] possible options: ["default", "genotype", "csq"]

    option_csq : options when selecting CSQ columns option="potential" or "important"

    options_genotype : when selecting genotype columns options=["potential"/"important", "all"/"common"]
    """

    data_path = Path(data_folder)

    assert len(samples) > 0, "Number of samples must be higher than 0"

    df = get_sample_data(str(data_path / samples[0]), file_types)
    if len(samples) > 1:
        for sample in samples[1:]:
            df2 = get_sample_data(str(data_path / sample), file_types)
            df = pd.concat([df, df2], ignore_index=True, axis=0)

    # Do default columns cleaning here
    if "default" in file_types:
        df = clean_default(df)
        
    # Do genotype columns cleaning here
    if "genotype" in file_types:
        df = drop_genotype_columns(df, options_genotype)
        
    # Do csq columns cleaning here
    if "csq" in file_types and len(file_types) > 1:
        df = drop_csq_columns(df, option_csq, csq_cols=True)
    elif "csq" in file_types:
        df = drop_csq_columns(df, option_csq)

    return df

def get_sample_data(sample_folder_path, file_types=None):
    """
    Returns raw csv converted data for a sample

    Usage example:

    EE_015 = get_sample_data("data/EE_015/", ["default", "csq"])
    
    EE_050 = get_sample_data("data/EE_050/", ["default", "csq"])
    
    EE_069 = get_sample_data("data/EE_069/", ["default", "csq"])
    
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
        return pd.read_csv(folder_path / f"{sample}_{file_types[0]}.csv.gz", sep=";", compression="gzip", low_memory=False).drop("Unnamed: 0", axis=1)
    
    df = pd.read_csv(folder_path / f"{sample}_{file_types[0]}.csv.gz", sep=";", compression="gzip", low_memory=False).drop("Unnamed: 0", axis=1)
    if file_types[0] == "csq":
        df = df.add_suffix("_csq")

    for f_type in file_types[1:]:
        df2 = pd.read_csv(folder_path / f"{sample}_{f_type}.csv.gz", sep=";", compression="gzip", low_memory=False).drop("Unnamed: 0", axis=1)
        if f_type == "csq":
            df2 = df2.add_suffix("_csq")
        df = pd.concat([df, df2], axis=1, ignore_index=False)
        
    return df

def clean_default(EE_default):
    EE_default.drop("ID", axis=1, inplace=True)
    unique_filters = set(';'.join(EE_default['FILTER']).split(';'))
    unique_filters = [f for f in unique_filters if f != "PASS" and f != "FAIL"]
    filter_cols = ["FILTER_" + f for f in unique_filters]
    for f_col, f_val in zip(filter_cols, unique_filters):
        EE_default[f_col] = EE_default['FILTER'].apply(lambda x: 1 if f_val in x else 0)

    EE_default_clean = EE_default.drop("FILTER", axis=1)
    return EE_default_clean

def drop_csq_columns(EE_csq, option, csq_cols=False):
    """
    option="potential" or "important"
    """
    drop_columns = ["TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "REFSEQ_MATCH",
                "SOURCE", "REFSEQ_OFFSET", "GIVEN_REF", "USED_REF", "BAM_EDIT", "DOMAINS", "HGVS_OFFSET", "AF", "AFR_AF",
                "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "cDNA_position", "CDS_position", "Protein_position"]
    
    drop_columns_gnomad = [c for c in EE_csq if c.startswith("gnomAD") and (c!="gnomADe_AF" and c!="gnomADg_AF")]

    if csq_cols:
        drop_columns = [c + "_csq" for c in drop_columns]
        drop_columns_gnomad = [c for c in EE_csq if c.startswith("gnomAD") and (c!="gnomADe_AF_csq" and c!="gnomADg_AF_csq")]

    EE_potential_csq = EE_csq.drop(drop_columns, axis=1)
    EE_potential_csq.drop(drop_columns_gnomad, axis=1, inplace=True)
    
    if option == "potential":
        return EE_potential_csq

    elif option == "important":
        potential_drop_columns = ["Consequence", "IMPACT", "CANONICAL", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SIFT", "PolyPhen",
                          "CLIN_SIG", "EVE_CLASS", "EVE_SCORE", "CADD_PHRED", "CADD_RAW", "LOEUF", "NMD", "SpliceAI_pred_DP_AG",
                          "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG",
                          "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "SpliceAI_pred_SYMBOL"]

        if csq_cols:
            potential_drop_columns = [c + "_csq" for c in potential_drop_columns]

        EE_important_csq = EE_potential_csq.drop(potential_drop_columns, axis=1)
        return EE_important_csq

    else:
        raise Exception("Option chosen incorrectly please select option='potential' or 'important'")
    
def drop_genotype_columns(EE_genotype, options):
    """
    options=["potential"/"important", "all"/"common"]

    - potential - potentially interesting columns are included in the final dataset
    - important - only the most promising columns are included in the final dataset
    - common - only common columns between the EE_015/EE_050 and EE_069 datasets are included in the final dataset
    - all - all columns in the EE_015/EE_050 and EE_069 datasets are included in the final dataset. 
    Choose this option if 069 is not included in the dataset.
    """
    assert (options[0] == "potential" or options[0] == "important"), "Incorrect options chosen please take look at the function description"
    assert (options[1] == "all" or options[1] == "common"), "Incorrect options chosen please take look at the function description"

    excess_069_columns = ['AC', 'AF', 'AN', 'BaseQRankSum', 'ClippingRankSum', 'ExcessHet', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR']
    excess_069_columns = [c for c in excess_069_columns if c in EE_genotype.columns]
    EE_genotype.drop(excess_069_columns, axis=1, inplace=True)

    drop_columns_acmg_amp =  [c for c in EE_genotype.columns if c.startswith("ACMG") or c.startswith("AMP") and c!="ACMG_class"]
    drop_columns_gnomad = [c for c in EE_genotype.columns if c.startswith("gnomad") and c != "gnomadExomes_AF" and c!= "gnomadGenomes_AF"]

    EE_genotype.drop("DP", axis=1, inplace=True)
    if "MMQ" in EE_genotype.columns:
        EE_genotype.drop("MMQ", axis=1, inplace=True)  
    EE_genotype.drop(drop_columns_acmg_amp, axis=1, inplace=True)
    EE_genotype.drop(drop_columns_gnomad, axis=1, inplace=True)

    if options[0] == "important":
        potential_drop_columns = ["ClinVarClass", "ClinVarDisease", "DANN_score", "MutationTaster_pred", "MutationTaster_score", "SIFT_score"]
        EE_genotype.drop(potential_drop_columns, axis=1, inplace=True)

    if options[1] == "common":
        different_columns = ['AS_FilterStatus', 'AS_SB_TABLE', 'ECNT', 'GERMQ', 'MBQ', 'MFRL', 'MPOS', 'POPAF', 'RPA', 'RU', 'STR', 'STRQ', 'TLOD', 'cosmicFathMMPrediction', 'cosmicFathMMScore']
        EE_genotype.drop(different_columns, axis=1, inplace=True)

    return EE_genotype