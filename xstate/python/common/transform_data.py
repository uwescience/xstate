""" Functions that transform DataFrames. """

from common import constants as cn
from common.data_provider import DataProvider

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
  """
  if df is None:
    if provider is None:
      provider = DataProvider()
      provider.do()
    df = provider.df_normalized
  df_trinary = makeTrinaryData(df, is_include_nan=False)
  dfg = df_trinary.groupby(df_trinary.columns.tolist())
  groups = dfg.groups
  data = {}
  for key, genes in groups.items():
    label = "--".join(genes)
    data[label] = list(key)
  df = pd.DataFrame(data)
  df_result = df.T
  df_result.columns = df_trinary.columns
  return df_result

def trinaryReadsDF(csv_file=None, df_sample=None,
    csv_dir=cn.SAMPLES_DIR, is_display_errors=True):
  """
  Creates trinary values for read counts w.r.t. data provider.
  (a) adjusting for gene length, (b) library size,
  (c) log2, (d) ratio w.r.t. T0.
  Data may come from an existing dataframe or a CSV file.
  :param str csv_file: File in "samples" directory.
      columns are: "GENE_ID", instance ids
  :param pd.DataFrame df_sample: columns are genes,
      index are instances, values are raw readcounts
  :param str csv_dir: directory where csv file is found
  :return pd.DataFrame: columns are genes, 
      indexes are instances, trinary values
  At least one of df_sample and csv_file must be non-null
  """
  provider = DataProvider(is_display_errors=is_display_errors)
  provider.do()
  if df_sample is None:
    path = os.path.join(csv_dir, csv_file)
    df_sample = pd.read_csv(path)
    df_sample.index = df_sample['GENE_ID']
    del df_sample['GENE_ID']
  #
  df_normalized = provider.normalizeReadsDF(df_sample)
  # Compute trinary values relative to original reads
  df_ref = sum(provider.dfs_adjusted_read_count)  \
      / len(provider.dfs_adjusted_read_count)  # Mean values
  ser_ref = df_ref[cn.REF_TIME]
  return calcTrinaryComparison(df_normalized, ser_ref=ser_ref)

def calcTrinaryComparison(df, ser_ref=None, threshold=1, is_convert_log2=True):
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
  MINVAL = 1e-12
  if is_convert_log2:
    if not ser_ref is None:
      ser_ref_log = ser_ref.apply(lambda v: np.log2(v))
    df_log = df.applymap(lambda v: np.log2(v)
        if v > MINVAL else np.log2(MINVAL))
  else:
    ser_ref_log = ser_ref
    df_log = df
  #
  if ser_ref is None:
    ser_ref_log = pd.Series(np.repeat(0, len(df)), index=df.index)
  #
  df_comp = pd.DataFrame()
  for col in df.columns:
    df_comp[col] = df_log[col] - ser_ref_log
  df_result = df_comp.applymap(
      lambda v: 0 if np.abs(v) < threshold else -1 if v < 0 else 1)
  return df_result
  
