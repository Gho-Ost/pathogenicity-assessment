import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.feature_selection import mutual_info_classif
from sklearn.preprocessing import LabelEncoder
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import classification_report
from sklearn.utils.class_weight import compute_class_weight
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier, BaggingClassifier, VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import accuracy_score
import dalex as dx
from sklearn.model_selection import LeaveOneOut
from sklearn.naive_bayes import GaussianNB
from EDA_module import *


def drop_CSQ(dfs):
    for df in dfs:
        df.drop(columns="CSQ", inplace=True)

def drop_ACMG_AMP(df):
    columns_to_drop = [col for col in df.columns if 'ACMG' in col and col != 'ACMG_class' or 'AMP' in col]
    return df.drop(columns=columns_to_drop)

def save_NaNs(df, output_file_path):
    nan_columns = df.columns[df.isnull().any()].tolist()
    nan_counts = {column: df[column].isnull().sum() for column in nan_columns}
    with open(output_file_path, 'w') as file:
        for column, count in nan_counts.items():
            file.write(f'{column}: {count} NaN values\n')

    print(f'Information saved to {output_file_path}')

def plot_columns_above_Nan_threshold(df, nan_threshold):
    columns_with_nans = []
    nan_counts = []

    for column in df.columns:
        nan_count = df[column].isna().sum()
        if nan_count > nan_threshold:
            columns_with_nans.append(column)
            nan_counts.append(nan_count)

    total_rows = len(df)

    plt.figure(figsize=(40, 6))
    plt.bar(columns_with_nans, nan_counts, color='skyblue')
    plt.axhline(y=total_rows, color='r', linestyle='--', label='Number of Rows')
    plt.xlabel('Columns')
    plt.ylabel('Count of NaN Values')
    plt.title('Comparison of Columns with More than 100,000 NaNs to Total Number of Rows')
    plt.legend()
    plt.xticks(rotation=45)
    plt.show()

def drop_Nans_threshold(df, nan_threshold):
    return df.dropna(axis=1, thresh=len(df) - nan_threshold)

def show_unique_vals(df):
    for column in df.select_dtypes(include=['object']):
        unique_values_count = df[column].nunique()
        print(f"Column '{column}' has {unique_values_count} unique values.")

def change_obj_to_category(df, unique_value_threshold):
    for column in df.select_dtypes(include=['object']):
        unique_values_count = df[column].nunique()
        if unique_values_count < unique_value_threshold:
            df[column] = df[column].astype('category')

def change_to_category(df, unique_value_threshold):
    for column in df.columns:
        unique_values_count = df[column].nunique()
        if unique_values_count < unique_value_threshold:
            df[column] = df[column].astype('category')

def map_objects_cats_to_nums(df):
    category_columns = df.select_dtypes(include=['category'])
    for column in category_columns:
        category_to_numeric_mapping = {category: index for index, category in enumerate(df[column].cat.categories)}
        df[column] = df[column].cat.codes

    object_columns = df.select_dtypes(include=['object'])
    for column in object_columns:
        df[column] = df[column].astype('category').cat.codes

def get_ACMG_mapping(df):
    acmg_class_names = list(df["ACMG_class"])
    acmg_class_vals = list(df["ACMG_class"])
    acmg_mapping={}
    for i in range(len(acmg_class_names)):
        name = acmg_class_names[i]
        if name not in list(acmg_mapping.keys()):
            acmg_mapping[name] = acmg_class_vals[i]
        if len(list(acmg_mapping.keys())) == 5:
            break
    return acmg_mapping

def create_infromation_plot(df, imputer_startegy="median"):
    target_column = 'ACMG_class'
    df_copy = df.copy()
    non_numeric_columns = df_copy.select_dtypes(exclude=['number']).columns
    label_encoders = {}

    for column in non_numeric_columns:
        label_encoders[column] = LabelEncoder()
        df_copy[column] = label_encoders[column].fit_transform(df_copy[column])

    imputer = SimpleImputer(strategy=imputer_startegy)
    X = df_copy.drop(columns=[target_column])
    X_imputed = imputer.fit_transform(X)
    y = df_copy[target_column]
    mutual_info = mutual_info_classif(X_imputed, y)


    mutual_info_df = pd.DataFrame({'Feature': X.columns, 'Mutual Information': mutual_info})
    mutual_info_df = mutual_info_df.sort_values(by='Mutual Information', ascending=False)
    top_features = mutual_info_df.head(15)

    plt.figure(figsize=(10, 6))
    plt.barh(top_features['Feature'], top_features['Mutual Information'], color='skyblue')
    plt.xlabel('Mutual Information')
    plt.ylabel('Feature')
    plt.title('Mutual Information with ACMG_class for Top 15 Features')
    plt.show()
