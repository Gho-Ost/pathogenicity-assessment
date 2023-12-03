import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.feature_selection import mutual_info_classif
from sklearn.preprocessing import LabelEncoder
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split


def get_target_plot(df, target):
    '''
    Plot proportions of value present in a chosen column
    '''
    plt.figure(figsize=(10, 6))
    sns.countplot(data=df, x=target, order=df[target].value_counts().index)
    plt.xlabel(target)
    plt.ylabel('Count')
    plt.title(f'Count of Values in {target}')
    plt.xticks(rotation=90)
    plt.show()

def get_nans_plot(df, nan_threshold = 100000):
    '''
    Plot numbers of nans for all columns with nans above threshold
    '''
    columns_with_nans = []
    nan_counts = []

    for column in df.columns:
        nan_count = df[column].isna().sum()
        if nan_count > nan_threshold:
            columns_with_nans.append(column)
            nan_counts.append(nan_count)

    total_rows = len(df)

    plt.figure(figsize=(40, 6))
    plt.bar(columns_with_nans, nan_counts)
    plt.axhline(y=total_rows, color='r', linestyle='--', label='Number of Rows')
    plt.xlabel('Columns')
    plt.ylabel('Count of NaN Values')
    plt.title('Comparison of Columns with More than 100,000 NaNs to Total Number of Rows')
    plt.legend()
    plt.xticks(rotation=45)
    plt.show()

def inspect_columns(df):
    '''
    A helper function that does a better job than df.info() and df.describe()
    '''
    result = pd.DataFrame({
        'unique': df.nunique() == len(df),
        'cardinality': df.nunique(),
        'with_null': df.isna().any(),
        'null_count': df.isnull().sum(),
        'null_pct': round((df.isnull().sum() / len(df)) * 100, 2),
        '1st_row': df.iloc[0],
        'random_row': df.iloc[np.random.randint(low=0, high=len(df))],
        'last_row': df.iloc[-1],
        'dtype': df.dtypes
    })
    return result

def plot_corr_heatmap_nans(df):
    '''
    Plot correlation between occurences of Nan values in a dataframe.
    '''
    columns_with_nan = df.columns[df.isnull().any()]
    df_with_nan = df[columns_with_nan]
    nan_masks = df_with_nan.isnull()
    nan_masks = nan_masks.astype(int)
    nan_correlation = nan_masks.corr()
    nan_correlation = nan_correlation.round(2)
    plt.figure(figsize=(40, 30))
    sns.heatmap(nan_correlation, annot=True, cmap='coolwarm', linewidths=.5)
    plt.title('Correlation Nan occurrences')
    plt.show()

def plot_corr_heatmap_nans_target(df, target):
    '''
    Plot correlation between occurences of Nan values in a dataframe including not-Nan target.
    '''
    columns_with_nan = df.columns[df.isnull().any()]
    columns_to_include = list(columns_with_nan) + [target]

    df_with_nan = df[columns_to_include]
    df_with_nan[target] = df_with_nan[target].astype('category').cat.codes
    nan_masks = df_with_nan.isnull()
    nan_masks = nan_masks.astype(int)
    nan_masks[target] = df_with_nan[target]
    nan_correlation = nan_masks.corr()
    nan_correlation = nan_correlation.round(2)
    plt.figure(figsize=(40, 30))
    sns.heatmap(nan_correlation, annot=True, cmap='coolwarm', linewidths=.5)
    plt.title('Correlation Nan occurrences')
    plt.show()

def get_nan_corr_treshold(df, target, correlation_threshold = 0.15):
    '''
    Get Columns with Nan correlation with target above a threshold
    '''
    columns_with_nan = df.columns[df.isnull().any()]
    columns_to_include = list(columns_with_nan) + [target]

    df_with_nan = df[columns_to_include]
    df_with_nan[target] = df_with_nan[target].astype('category').cat.codes
    nan_masks = df_with_nan.isnull()
    nan_masks = nan_masks.astype(int)
    nan_masks[target] = df_with_nan[target]
    nan_correlation = nan_masks.corr()
    nan_correlation = nan_correlation.round(2)
    correlated_columns = nan_correlation[target][abs(nan_correlation[target]) > correlation_threshold].index
    return correlated_columns

def plot_dtypes(df):
    '''
    Plot proportions of data types present in the dataframe
    '''
    dtype_counts={}
    for val in df.dtypes:
        try:
            dtype_counts[str(val)]+=1
        except:
            dtype_counts[str(val)]=1

    data_types = list(dtype_counts.keys())
    counts = list(dtype_counts.values())

    plt.figure(figsize=(8, 6))
    plt.bar(data_types, counts)
    plt.xlabel('Data Types')
    plt.ylabel('Count')
    plt.title('Counts of Different Data Types')
    plt.xticks(rotation=45)
    plt.show()

def preprocess(df, one_hot_nans = True, fill_median = False):
    '''
    Preprocess the datframe according to the proposed schema
    '''
    new_df = df.copy()
    if one_hot_nans:
        for column in new_df.columns:
            if new_df[column].isnull().any():
                new_df[column + '_is_nan'] = new_df[column].isnull().astype(int)
    for column in new_df.columns:
        if new_df[column].isnull().any():
            if new_df[column].dtype == 'object':
                new_df[column].fillna("No Value", inplace=True)
            else:
                if fill_median:
                    new_df[column].fillna(new_df[column].median(), inplace=True)
                else:
                    new_df[column].fillna(new_df[column].min() - 1, inplace=True)
    return new_df

def encode(df, target_column, target_mapping):
    '''
    Encode the datframe according to the proposed schema
    '''
    new_df = df.copy()
    label_encoders = {}
    for column in new_df.columns:
        if new_df[column].dtype == 'object' and column!=target_column:
            label_encoders[column] = LabelEncoder()
            try:
                new_df[column] = label_encoders[column].fit_transform(new_df[column])
            except:
                new_df[column] = label_encoders[column].fit_transform(new_df[column].astype(str))
        elif new_df[column].dtype == 'object' and column==target_column:
            new_df[target_column] = new_df[target_column].map(target_mapping)
    return new_df, label_encoders, target_mapping

def get_mapping(encoders, expected_column):
    '''
    Get mapping from encoders after performing the encode() function
    '''
    for column, encoder in encoders.items():
        if column == expected_column:
            print(f"{column} mapping: {dict(zip(encoder.classes_, encoder.transform(encoder.classes_)))}")
            break

def get_corr_heatmap(df, k=20, plot_size = (20, 10)):
    '''
    Plot correlation heatmap
    '''
    numeric_df = df.select_dtypes(include=['int64', 'float64', 'int32'])

    correlation_matrix = numeric_df.corr()
    average_correlation = correlation_matrix.abs().mean(axis=1)
    average_correlation = average_correlation.round(2)
    top_features = average_correlation.nlargest(k).index

    plt.figure(figsize=plot_size)
    sns.heatmap(correlation_matrix.loc[top_features, top_features], annot=True, cmap='coolwarm', linewidths=.5)
    plt.title(f'Correlation Heatmap for Top {k} Features with Highest Average Correlation')
    plt.show()

def get_corr_target(df, target_column, k=20, plot_size = (10, 6)):
    '''
    Get correlation with target and plot top k atributes with highest correlation
    '''
    correlation_with_ACMG_class = df.corr()[target_column].drop(target_column)

    top_columns = correlation_with_ACMG_class.abs().nlargest(k).index

    plt.figure(figsize=plot_size)
    # colors = ['#67a9cf' if c >= 0 else '#ef8a62' for c in correlation_with_ACMG_class[top_columns]]
    colors = ['#74add1' if c >= 0 else '#f46d43' for c in correlation_with_ACMG_class[top_columns]]
    
    sns.barplot(x=correlation_with_ACMG_class[top_columns], y=top_columns, orient='h', palette=colors)
    plt.xlabel(f'Correlation with {target_column}')
    plt.ylabel('Column Name')
    plt.title(f'Correlation of Top {k} Columns with {target_column}')
    plt.show()

def get_mutual_info_plot(df, target_column, k=20, plot_size = (10, 6)):
    '''
    Get mutual information with target and plot top k atributes with highest mutual information
    '''
    df_copy = df.copy()
    X = df_copy.drop(columns=[target_column])
    y = df_copy[target_column]
    mutual_info = mutual_info_classif(X, y)


    mutual_info_df = pd.DataFrame({'Feature': X.columns, 'Mutual Information': mutual_info})
    mutual_info_df = mutual_info_df.sort_values(by='Mutual Information', ascending=False)
    top_features = mutual_info_df.head(k)

    plt.figure(figsize=plot_size)
    plt.barh(top_features['Feature'], top_features['Mutual Information'])
    plt.xlabel('Mutual Information')
    plt.ylabel('Feature')
    plt.title(f'Mutual Information with {target_column} for Top {k} Features')
    plt.show()