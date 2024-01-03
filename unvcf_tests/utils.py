"""
Module for custom utility functions.
"""

import pandas as pd
import numpy as np

from pathlib import Path
from warnings import warn

def get_dataset(data_folder, samples, file_type, option_csq=None, options_genotype=None, with_default=False):
    """
    Returns cleaned up dataset from chosen samples and file types. 
    Make sure that a correct folder structure is prepared before using this function.

    Parameters
    ----------
    data_folder : path to data folder e.g. "./data"

    samples : samples to be included in the dataset e.g. ["EE_015", "EE_050"] possible options: ["EE_015", "EE_050", "EE_069"]

    file_type : file type to be included in the dataset possible options: "genotype"/"csq"/"both"

    option_csq : options when selecting CSQ columns possible options: "potential"/"important"

    options_genotype : when selecting genotype columns possible options : ["potential"/"important", "all"/"common"]

    with_default : include "default" columns
    """

    data_path = Path(data_folder)

    assert len(samples) > 0, "Number of samples must be higher than 0"
    
    if file_type == "csq" or file_type == "both":
        assert option_csq != None, "CSQ options cannot be None"

    if file_type == "genotype" or file_type == "both":
        assert options_genotype != None, "Genotype options cannot be None"

    # Load default info for all chosen samples
    if with_default:
        default_dataframes = []

        for sample in samples:
            df = pd.read_csv(data_path / sample / f"{sample}_default.csv.gz", sep=";", compression="gzip", low_memory=False).drop("Unnamed: 0", axis=1)
            default_dataframes.append(df)
    
        df_default = pd.concat(default_dataframes, ignore_index=True, axis=0)
        
        # Clean up default info
        df_default_clean = clean_default_columns(df_default)

    # Load genotype for all chosen samples
    genotype_dataframes = []
    for sample in samples:
        df = pd.read_csv(data_path / sample / f"{sample}_genotype.csv.gz", sep=";", compression="gzip", low_memory=False).drop("Unnamed: 0", axis=1)
        genotype_dataframes.append(df)

    df_genotype = pd.concat(genotype_dataframes, ignore_index=True, axis=0)

    # If genotype data is required, clean it up
    if file_type == "genotype" or file_type == "both":
        df_genotype_filtered = drop_genotype_columns(df_genotype, options_genotype)
        df_genotype_clean = clean_genotype_columns(df_genotype_filtered)

    # Load csq for all chosen samples
    if file_type == "csq" or file_type == "both":
        csq_dataframes = []

        for sample in samples:
            df = pd.read_csv(data_path / sample / f"{sample}_csq.csv.gz", sep=";", compression="gzip", low_memory=False).drop("Unnamed: 0", axis=1)
            csq_dataframes.append(df)

        df_csq = pd.concat(csq_dataframes, ignore_index=True, axis=0)

        # Clean up csq data
        df_csq_filtered = drop_csq_columns(df_csq, option_csq)
        df_csq_clean = clean_csq_columns(df_csq_filtered)

        # Add target to csq data
        if file_type == "csq":
            df_csq_clean["ACMG_class"] = df_genotype["ACMG_class"]

    # Select proper dataframe to be returned
    if file_type == "both":
        df = pd.concat([df_genotype_clean, df_csq_clean], axis=1, ignore_index=False)

    elif file_type == "csq":
        df = df_csq_clean

    elif file_type == "genotype":
        df = df_genotype_clean

    # Add default data
    if with_default:
        df = pd.concat([df_default_clean, df], axis=1, ignore_index=False)

    return df

def clean_default_columns(df_default):
    df_default.drop("ID", axis=1, inplace=True)
    unique_filters = set(';'.join(df_default['FILTER']).split(';'))
    unique_filters = [f for f in unique_filters if f != "PASS" and f != "FAIL"]
    filter_cols = ["FILTER_" + f for f in unique_filters]
    for f_col, f_val in zip(filter_cols, unique_filters):
        df_default[f_col] = df_default['FILTER'].apply(lambda x: 1 if f_val in x else 0)

    df_default_clean = df_default.drop("FILTER", axis=1)

    # Add relative POS
    chromosome_lengths = [0, 16569, 248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 
                      145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 
                      90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895]
    summed_lengths = np.cumsum(chromosome_lengths, dtype=np.int64)

    chromosome_names = [f"chr{i}" for i in range(0, 23)]
    chromosome_names[0] = "chrM"
    chromosome_names.append("chrX")
    chromosome_names.append("chrY")

    chromosome_lengths_dict = dict([n, l] for n, l in zip(chromosome_names, summed_lengths))
    df_default_clean['POS_REL'] = df_default_clean.apply(lambda row: row['POS'] + chromosome_lengths_dict.get(row['#CHROM'], 0), axis=1)

    # NaNs in unusual chromosomes POS_REL
    df_default_clean.loc[~df_default_clean['#CHROM'].isin(chromosome_names), 'POS_REL'] = np.nan

    return df_default_clean

def _extract_MutationTaster_values(row):
    # MutationTaster_pred and MutationTaster_score split select in order {A, D, P, N}>0.5 and with highest confidence otherwise highest confidence
    if row['MutationTaster_pred'] == '0':
        return pd.Series([np.nan, np.nan], index=["MutationTaster_selected_pred", "MutationTaster_selected_score"])

    symbols = ['A', 'D', 'P', 'N']
    pred_symbols = row['MutationTaster_pred'].split('%3B')
    score_values = [float(x) for x in row['MutationTaster_score'].split('%3B')]

    top_score = -1
    top_symbol = None

    # Check for symbols A, D, P, N with score above 0.5 and select the one with the highest score
    for symbol in symbols:
        if symbol in pred_symbols:
            index_symbol = pred_symbols.index(symbol)
            selected_symbol = symbol
            selected_score = score_values[index_symbol]
            
            if selected_score >= 0.5:
                return pd.Series([selected_symbol, selected_score], index=["MutationTaster_selected_pred", "MutationTaster_selected_score"])
            
            elif selected_score > top_score:
                top_score = selected_score
                top_symbol = selected_symbol

    return pd.Series([top_symbol, top_score], index=["MutationTaster_selected_pred", "MutationTaster_selected_score"])
    

def clean_genotype_columns(df_genotype):
    """
    Fixes genotype columns with array values.
    """
    object_columns = set(df_genotype.select_dtypes(object).columns)
    non_arrays = ["RU", "cosmicFathMMPrediction", "ACMG_class"]

    # Take undropped array columns
    separable_columns = [c for c in object_columns if c not in non_arrays and c in df_genotype.columns]

    # Clean up columns
    for c in separable_columns:
        # MutationTaster_pred and _score
        if c == "MutationTaster_pred":
            df_genotype[["MutationTaster_selected_pred", "MutationTaster_selected_score"]] = df_genotype.apply(_extract_MutationTaster_values, axis=1)
            df_genotype.drop("MutationTaster_pred", axis=1, inplace=True)
            df_genotype.drop("MutationTaster_score", axis=1, inplace=True)
        
        # CGDinheritance split - if '0' then 0 else 1
        elif c == "CGDinheritance":
            df_genotype["CGDinheritance_exist"] = df_genotype["CGDinheritance"].apply(lambda x: 1 if x!='0' else 0)
            df_genotype.drop("CGDinheritance", axis=1, inplace=True)

        # function split - take possible results and one hot encode
        elif c == "function":
            functions = ["0", "NMD", "3'utr", "5'utr", "3'flank", "5'flank", "coding", "non-coding%40exon", "intronic", "splicing", "splicing-ACMG"]
            function_cols = ["function_" + c for c in functions]

            for function, col in zip(functions, function_cols):
                df_genotype[col] = df_genotype["function"].apply(lambda x: 1 if function in str(x) else 0)

            df_genotype.drop("function", axis=1, inplace=True)

        # coding_impact split - one hot encode for possible values
        elif c == "coding_impact":
            coding_impact_unique = set([])
            for c in df_genotype["coding_impact"].unique():
                if not isinstance(c, float):        
                    for _c in c.split(","):
                        coding_impact_unique.add(_c)
                else:
                    coding_impact_unique.add(c)

            coding_impact_cols = ["coding_impact_" + c for c in coding_impact_unique]

            for element, col in zip(coding_impact_unique, coding_impact_cols):
                df_genotype[col] = df_genotype["coding_impact"].apply(lambda x: 1 if element in str(x) else 0)

            df_genotype.drop("coding_impact", axis=1, inplace=True)

    return df_genotype

def drop_genotype_columns(df_genotype, options):
    """
    Drops unnecessary genotype columns.
    options=["potential"/"important", "all"/"common"]

    - potential - potentially interesting columns are included in the final dataset
    - important - only the most promising columns are included in the final dataset
    - common - only common columns between the EE_015/EE_050 and EE_069 datasets are included in the final dataset
    - all - all columns in the EE_015/EE_050 and EE_069 datasets are included in the final dataset. 
    Choose this option if 069 is not included in the dataset.
    """
    assert (options[0] in ["potential", "important"]), "Incorrect options chosen please take look at the function description"
    assert (options[1] in ["all", "common"]), "Incorrect options chosen please take look at the function description"

    excess_069_columns = ['AC', 'AF', 'AN', 'BaseQRankSum', 'ClippingRankSum', 'ExcessHet', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR']
    excess_069_columns = [c for c in excess_069_columns if c in df_genotype.columns]
    df_genotype.drop(excess_069_columns, axis=1, inplace=True)

    drop_columns_acmg =  [c for c in df_genotype.columns if c.startswith("ACMG") and c!="ACMG_class"]
    drop_columns_amp =  [c for c in df_genotype.columns if c.startswith("AMP")]
    drop_columns_gnomad = [c for c in df_genotype.columns if c.startswith("gnomad") and c != "gnomadExomes_AF" and c!= "gnomadGenomes_AF"]

    # Drop these columns for both datasets
    df_genotype.drop(drop_columns_acmg, axis=1, inplace=True)   
    df_genotype.drop(drop_columns_gnomad, axis=1, inplace=True)

    df_genotype.drop("DP", axis=1, inplace=True)
    df_genotype.drop("Gene", axis=1, inplace=True)
    df_genotype.drop("hgvs", axis=1, inplace=True)
    df_genotype.drop("ClinVarClass", axis=1, inplace=True)

    # Drop columns from 015 and 050 datasets
    if "MMQ" in df_genotype.columns:
        df_genotype.drop("MMQ", axis=1, inplace=True)  
        df_genotype.drop("AS_SB_TABLE", axis=1, inplace=True)
        df_genotype.drop("AS_FilterStatus", axis=1, inplace=True)
        df_genotype.drop(drop_columns_amp, axis=1, inplace=True)   

    # If only important columns should be left
    if options[0] == "important":
        potential_drop_columns = ["ClinVarDisease", "DANN_score", "MutationTaster_pred", "MutationTaster_score", "SIFT_score"]
        
        # Make sure the columns still exist in the dataset
        potential_drop_columns = [c for c in potential_drop_columns if c in df_genotype.columns]
        df_genotype.drop(potential_drop_columns, axis=1, inplace=True)

    # If common columns of 015/050 and 069 should be kept
    if options[1] == "common":
        different_columns = ['ECNT', 'GERMQ', 'MBQ', 'MFRL', 'MPOS', 'POPAF', 'RPA', 'RU', 'STR', 'STRQ', 'TLOD', 'cosmicFathMMPrediction', 'cosmicFathMMScore']
        df_genotype.drop(different_columns, axis=1, inplace=True)

    return df_genotype

def clean_csq_columns(df_csq):
    """
    Fixes csq columns with array values.
    """
    # Add other columns here when a process to clean them up is known
    separable_columns = ["PHENOTYPES", "CLIN_SIG", "HGNC_ID", "Existing_variation", "FLAGS", "PUBMED", "SIFT", "PolyPhen"]

    # Take undropped columns
    separable_columns = [c for c in separable_columns if c in df_csq.columns]

    # Clean up columns
    for c in separable_columns:
        # PHENOTYPES split - test if a value exist in the cell
        if c == "PHENOTYPES" or c == "PUBMED" or c == "CLIN_SIG" or c == "HGNC_ID":
            df_csq[f"{c}_exist"] = df_csq[c].apply(lambda x: 1 if pd.notna(x) else 0)
            df_csq.drop(c, axis=1, inplace=True)
        
        # Existing_variation - split on: exists in COSV/rs
        elif c == "Existing_variation":
            existing_var_unique = set([])
            for c in df_csq["Existing_variation"].unique():
                if not isinstance(c, float):
                    for _c in c.split("&"):
                        existing_var_unique.add("".join([x for x in _c if not x.isnumeric()]))

            existing_var_cols = ["Existing_variation_" + x for x in list(existing_var_unique)]

            for prefix, col in zip(existing_var_unique, existing_var_cols):
                df_csq[col] = df_csq["Existing_variation"].apply(lambda x: 1 if prefix in str(x) else 0)

            df_csq.drop("Existing_variation", axis=1, inplace=True)

        # FLAGS - separate to start and end flag existance in row
        elif c == "FLAGS":
            flags = ['cds_start_NF', 'cds_end_NF']

            df_csq["FLAGS_start"] = df_csq["FLAGS"].apply(lambda x: 1 if flags[0] in str(x) else 0)
            df_csq["FLAGS_start"] = df_csq["FLAGS"].apply(lambda x: 1 if flags[1] in str(x) else 0)
            df_csq.drop("FLAGS", axis=1, inplace=True)
        
        # SIFT and PolyPhen splits - separate to a class and prediction
        elif c == "SIFT" or c == "PolyPhen":
            df_csq[f'{c}_class'] = df_csq[c].str.extract(r'([^\(]+)')
            df_csq[f'{c}_pred'] = df_csq[c].str.extract(r'\(([^)]+)\)').astype(float)
            df_csq.drop(c, axis=1, inplace=True)

    return df_csq

def drop_csq_columns(df_csq, option):
    """
    Drops unnecessary csq columns.
    option="potential" or "important"
    """
    assert (option in ["potential", "important"]), "Incorrect option selected please choose 'important' or 'potential'"

    # Drop the unnecessary columns
    drop_columns = ["TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "REFSEQ_MATCH", "Gene",
                "SOURCE", "REFSEQ_OFFSET", "GIVEN_REF", "USED_REF", "BAM_EDIT", "DOMAINS", "HGVS_OFFSET", "AF", "AFR_AF",
                "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "cDNA_position", "CDS_position", "Protein_position", "HGVSp", "HGVSc", "Consequence",
                "TRANSCRIPTION_FACTORS", "SOMATIC", "PHENO", "VAR_SYNONYMS"]
    
    drop_columns_gnomad = [c for c in df_csq if c.startswith("gnomAD") and (c!="gnomADe_AF" and c!="gnomADg_AF")]

    df_csq.drop(drop_columns, axis=1, inplace=True)
    df_csq.drop(drop_columns_gnomad, axis=1, inplace=True)

    # If only important columns should be evaluated
    if option == "important":
        potential_drop_columns = ["IMPACT", "CANONICAL", "SYMBOL", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SIFT", "PolyPhen",
                          "CLIN_SIG", "EVE_CLASS", "EVE_SCORE", "CADD_PHRED", "CADD_RAW", "LOEUF", "NMD", "SpliceAI_pred_DP_AG",
                          "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG",
                          "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "SpliceAI_pred_SYMBOL", "FLAGS"]

        # Make sure the columns are not yet dropped
        potential_drop_columns = [c for c in potential_drop_columns if c in df_csq.columns]
        
        df_csq.drop(potential_drop_columns, axis=1, inplace=True)
    
    return df_csq

def get_sample_data(sample_folder_path, file_types=None):
    """
    Returns raw csv converted data for a sample

    Usage example:

    EE_015 = get_sample_data("data/EE_015/", ["default", "csq"])
    
    EE_050 = get_sample_data("data/EE_050/", ["default", "csq"])
    
    EE_069 = get_sample_data("data/EE_069/", ["default", "csq"])
    
    df = pd.concat([EE_015, EE_050, EE_069], ignore_index=True, axis=0)
    """

    warn("This function is depracated watch out for changes in columns")

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