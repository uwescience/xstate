""" Functions that transform DataFrames. """

from common import constants as cn
from common.data_provider import DataProvider
from common import data_provider

import copy
import os
import pandas as pd
import numpy as np


# TODO: Should nan values be a trinary 0?
def makeTrinaryData(df=None, min_abs=1.0, is_include_nan=True):
  """
  Thresholds data based on its absolute magnitude.
  Values are assigned as -1, 0, 1
  :param pd.DataFrame df: default is provider.df_normalized
    values are in log2 units
  :param float min_abs: minimal absolute value to threshold.
  :param bool is_include_nan: Include nan values; else set to 0
  :return pd.DataFrame: same index and columns as df
  """
  if df is None:
    provider = DataProvider()
    provider.do()
    df = provider.df_normalized
  df_result = df.copy()
  df_result = df_result.applymap(
      lambda v: 0 if np.abs(v) < min_abs else -1 if v < 0 else 1)
  if is_include_nan:
      df_result = df_result.applymap(
         lambda v: np.nan if v==0 else v)
  return df_result

def aggregateGenes(df=None, provider=None):
  """
  Combines genes that are perfectly correlated in time for trinary
  values.
  :param DataFrame df: dataframe to transform
  :param DataProvider provider: uses df_normalized
  :return pd.DataFrame: names are combined for aggregated
      genes; calculates trinary values
      indexes are unchanged
  """
  if df is None:
    if provider is None:
      provider = DataProvider()
      provider.do()
    df = provider.df_normalized
  df_trinary = makeTrinaryData(df, is_include_nan=False)
  dfg = df_trinary.T.groupby(
      df_trinary.T.columns.tolist())
  groups = dfg.groups
  data = {}
  for key, genes in groups.items():
    label = cn.GENE_SEPARATOR.join(genes.values.tolist())
    data[label] = list(key)
  df = pd.DataFrame(data)
  df.index = df_trinary.index
  return df

def stripReplicaString(names):
  """
  Strips the replica information from a time.
  :param list-str names:
  :return list-str:
  """
  new_names = []
  for name in names:
    splits = name.split(data_provider.SEPARATOR)
    new_names.append(splits[0])
  return new_names

def readGeneCSV(csv_file, csv_dir=cn.SAMPLES_DIR):
  """
  Reads a CSV file with the column cn.GENE_ID
  :param str csv_file: File in "samples" directory.
      columns are: "GENE_ID", instance ids
  :param str csv_dir: directory where csv file is found
  :return pd.DataFrame:
      index: cn.GENE_ID
  """
  path = os.path.join(csv_dir, csv_file)
  df = pd.read_csv(path)
  df.index = df[cn.GENE_ID]
  del df[cn.GENE_ID]
  return df

def trinaryReadsDF(csv_file=None, df_sample=None,
    csv_dir=cn.SAMPLES_DIR, is_display_errors=True,
    ser_ref=None,
    is_normalized=False,
    is_time_columns=True, col_ref=None, is_convert_log2=True):
  """
  Creates trinary values for read counts w.r.t. data provider.
  (a) adjusting for gene length, (b) library size,
  (c) log2, (d) ratio w.r.t. T0. The T0 values is
  the average of the values in default DataProvider.
  Data may come from an existing dataframe or a CSV file.
  :param str csv_file: File in "samples" directory.
      columns are: "GENE_ID", instance ids
  :param pd.DataFrame df_sample: columns are genes,
      index are instances, values are raw readcounts
  :param pd.Series ser_ref: Reference values for
      calculating expression levels
  :param bool is_time_columns: a time column is present
  :param bool is_normalized: data are already normalized
  :param str csv_dir: directory where csv file is found
  :param str col_ref: column to use as reference in
      normalization
  :return pd.DataFrame: columns are genes, 
      indexes are instances, trinary values
  Exactly one of df_sample and csv_file must be non-null
  """
  provider = DataProvider(is_display_errors=is_display_errors)
  provider.do()
  # Get the sample data to transform 
  if df_sample is None:
    df_sample = readGeneCSV(csv_file, csv_dir=csv_dir)
  # Normalize the samples
  if is_normalized:
    df_normalized = df_sample
  else:
    df_normalized = provider.normalizeReadsDF(df_sample,
        is_time_columns=is_time_columns)
  # Construct the reference data
  if ser_ref is None:
    # Compute trinary values relative to original reads
    if col_ref is None:
      dfs = copy.deepcopy(provider.dfs_adjusted_read_count)
      for df in dfs:
        df.columns = stripReplicaString(df.columns)
      df_ref = sum(dfs) / len(provider.dfs_adjusted_read_count)
      col_name = provider.getT0s(df_ref)[0]
      ser_ref = df_ref[col_name]
    else:
      ser_ref = df_normalized[col_ref]
      del df_normalized[col_ref]
    df = calcTrinaryComparison(df_normalized, ser_ref=ser_ref,
        is_convert_log2=is_convert_log2)
  return df

def convertToLog2(obj):
  """
  Converts to log2 units.

  Parameters
  ----------
  obj: DataFrame / Series
  
  Returns
  -------
  DataFrame
  """
  if isinstance(obj, pd.DataFrame):
    return obj.applymap(lambda v: np.log2(v)
        if v > cn.MIN_VALUE else np.log2(cn.MIN_VALUE))
  else:  # Series
    return obj.apply(lambda v: np.log2(v)
        if v > cn.MIN_VALUE else np.log2(cn.MIN_VALUE))

def calcTrinaryComparison(df, ser_ref=None,
    threshold=1, is_convert_log2=True):
  """
  Calculates trinary values of a DataFrame w.r.t. a reference in
  log2 units.
  :param pd.Series ser_ref: reference values
  :param pd.DataFrame df: comparison values; columns are instances,
      has same inde as ser_ref
  :param float threshold: comparison threshold.
  :param bool is_convert_log2: convert to log2
  :return pd.DataFrame: trinary values resulting from comparisons
    -1: df is less than 2**threshold*ser_ref
     1: df is greater than 2**threshol*ser_ref
     0: otherwise
  """
  if is_convert_log2:
    if ser_ref is not None:
      ser_ref_log = convertToLog2(ser_ref)
    df_log = convertToLog2(df)
  else:
    ser_ref_log = ser_ref
    df_log = df
  #
  if ser_ref is None:
    ser_ref_log = pd.Series(np.repeat(0, len(df)), index=df.index)
  #
  df_comp = df_log.copy()
  # Find the common indices
  indices = set(df_log.index).intersection(ser_ref_log.index)
  df_log = df_log.loc[indices, :]
  ser_ref_log = ser_ref_log[indices]
  df_comp_T = df_log.T - ser_ref_log
  # Drop the nan columns, those genes for which there is no reference
  df_comp = (df_comp_T.dropna(axis=1, how='all')).T
  df_result = df_comp.applymap(
      lambda v: 0 if np.abs(v) < threshold else -1 if v < 0 else 1)
  return df_result

def removeGenesWithExcessiveReplicationVariance(df_X, max_var=None):
  """
  Removes Genes with excessive variation in the variance of their
  trinary values. Assumes a single digit replication.

  Parameters
  ----------
  df_X: DataFrame
      index: str <instance>.replication-digit
         ex: T10.2
  max_var: float
  
  Returns
  -------
  DataFrame (Trinary features)
  """
  df = df_X.copy()
  if max_var is None:
    return df
  #
  df.index = [i[0:-2] for i in df_X.index]
  df = df.sort_index()
  ser = df.groupby(df.index).std().sum()
  ser = ser.sort_values()
  ser_sub = ser[ser <= max_var]
  columns = list(ser_sub.index)
  return df_X[columns]
