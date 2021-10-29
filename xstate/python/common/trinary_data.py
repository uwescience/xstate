"""Normalized and Trinary Data."""

"""
df_X - feature matrix values with numeric column names.
       T0 is deleted
ser_y - state values
features - names of genes
state_dct - dictionary mapping string names to numbers
"""

import common.constants as cn
from common.data_provider import DataProvider
from common import data_provider
import common.transform_data as transform_data
from common_python.classifier import util_classifier

import collections
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

T1_INDEX = "T1"
MIN_NUM_NORMOXIA = 2  # Minimum number of normoxia states
FILE_GALAGAN = "galagan_raw_hypoxia_ts.csv" 
FILE_AM_MDM = "AM_MDM_Mtb_transcripts_DEseq.csv"
FILE_AW = "AW_plus_v_AW_neg_Mtb_transcripts_DEseq.csv"
FILE_RUSTAD = "rustad_hypoxia_dataset_GSE9331.csv"
FILE_GSE167232 = "GSE167232_mtb_transcriptome_counts_normalized_filtered.csv"
SHERMAN_INDUCED_PATH = os.path.join(cn.SAMPLES_DIR,
    "sherman_induced_mtb.txt")
SHERMAN_REPRESSED_PATH = os.path.join(cn.SAMPLES_DIR,
    "sherman_repressed_mtb.txt")
SAMPLES = ["AM_MDM", "AW", "sherman", "galagan", "rustad", "GSE167232"]
SampleData = collections.namedtuple("SampleData",
    SAMPLES)


################## FUNCTIONS ###############
def _subsetToRegulators(df_X):
  provider = DataProvider()
  provider.do()
  regulators = provider.df_trn_unsigned[cn.TF]
  regulators = list(set(regulators))
  keys = set(df_X.columns).intersection(
      regulators)
  return df_X[list(keys)]

def _getTrinaryFromGeneLists(
    repressed_path=SHERMAN_REPRESSED_PATH,
    induced_path=SHERMAN_INDUCED_PATH):
  """
  Creates a feature vector from a list of indiced
  and repressed genes.

  Parameters
  ----------

  repressed_path: str
  induced_path: str

  Return
  ______

  pd.DataFrame with single column
      index: gene
      value: trinary     
  """
  provider = DataProvider()
  provider.do()
  genes = provider.dfs_read_count[0].index.tolist()
  #
  def get_list(path):
    if path is None:
      return []
    with open(path, "r") as fd:
      items = [l.strip() for l in fd.readlines()]
    return list(set(genes).intersection(items))
    #
  represseds = get_list(repressed_path)
  induceds = get_list(induced_path)
  ser = pd.Series(np.repeat(0, len(genes)))
  ser.index = genes
  ser.loc[represseds] = -1
  ser.loc[induceds] = 1
  return pd.DataFrame(ser)

def getSampleData(is_regulator=True,
    is_display_errors=False, is_bioreactor_ref=True):
  """
  Acquires data obtain from other soruces.
  :param bool is_regulator: use regulators for TRN
  :param bool is_bioreactor_ref: use the DataProvider
      T0 data as the reference data
  :return SampleData. Each element is pd.DataFrame:
      columns: feature
      index: condition
      values: trinary
  """
  def makeSamples(csv_file, ref_sel_func, is_time_columns, df_data=None):
    """
    Constructs samples with respect to reference instances within the same
    data set.

    Parameters
    ----------
    csv_file: str
        CSV file for the sample data
    ref_sel_func: function
        parameters: dataframe index (str)
        returns: bool (include if True)
    is_time_columns: bool
        has a "time" column
    df_data: dataframe
        count to be converted into TrinaryData
    Returns
    -------
    DataFrame: Trinary values of samples
    """
    if df_data is None:
      df_data = transform_data.readGeneCSV(csv_dir=cn.SAMPLES_DIR,
          csv_file=csv_file).T
    provider = DataProvider(is_display_errors=False)
    provider.do()
    df_normalized = provider.normalizeReadsDF(
        df_data.T, is_time_columns=False).T
    ref_idxs = [i for i in df_normalized.index if ref_sel_func(i)]
    df_ref = df_normalized.loc[ref_idxs, :]
    ser_ref = df_ref.mean()
    sample_idxs = list(set(df_normalized.index).difference(ref_idxs))
    df_sample_raw = df_normalized.loc[sample_idxs, :]
    df_sample_trinary = transform_data.calcTrinaryComparison(df_sample_raw.T,
        ser_ref=ser_ref).T
    return df_sample_trinary
  #
  # AM/MDM
  if is_bioreactor_ref:
    df_AM_MDM = transform_data.trinaryReadsDF(
        is_display_errors=is_display_errors,
        csv_file=FILE_AM_MDM,
        is_time_columns=False).T
  else:
    ref_sel_func = lambda i: ("AM" in i) and (not "1" in i)
    df_AM_MDM = makeSamples(FILE_AM_MDM, ref_sel_func, False)
  if is_regulator:
    df_AM_MDM = _subsetToRegulators(df_AM_MDM)
  # AW
  if is_bioreactor_ref:
    df_AW = transform_data.trinaryReadsDF(
        csv_file=FILE_AW,
        is_display_errors=is_display_errors,
        is_time_columns=False).T
  else:
    ref_sel_func = lambda i: ("neg" in i) and (not "1" in i)
    df_AW = makeSamples(FILE_AW, ref_sel_func, False)
  if is_regulator:
    df_AW = _subsetToRegulators(df_AW)
  # Galagn data
  if is_bioreactor_ref:
    df_galagan = _getGalaganData(
        is_display_errors=is_display_errors, is_trinary=True)
  else:
    df_data = _getGalaganData(
        is_display_errors=is_display_errors, is_trinary=False)
    ref_sel_func = lambda i: ("d1." in i) and ("rep1" not in i)
    df_galagan = makeSamples(None, ref_sel_func, False, df_data=df_data)
  if is_regulator:
    df_galagan = _subsetToRegulators(df_galagan)
  df_sherman = _getTrinaryFromGeneLists()
  df_sherman = df_sherman.transpose()
  if is_regulator:
    df_sherman = _subsetToRegulators(df_sherman)
  # Add Rustad data
  if is_bioreactor_ref:
    df_rustad = transform_data.trinaryReadsDF(
        csv_file=FILE_RUSTAD,
        is_display_errors=is_display_errors,
        is_time_columns=False).T
  else:
    ref_sel_func = lambda i: ("_4hr_" in i) and ("rep6" not in i)
    df_rustad = makeSamples(FILE_RUSTAD, ref_sel_func, False)
  if is_regulator:
    df_rustad = _subsetToRegulators(df_rustad)
  # GSE167232
  #
  df_GSE167232 = transform_data.trinaryReadsDF(
      csv_file=FILE_GSE167232,
      is_display_errors=is_display_errors,
      is_normalized=True,
      is_time_columns=False).T
  if is_regulator:
    df_GSE167232 = _subsetToRegulators(df_GSE167232)
  #
  sample_data = SampleData(
      AM_MDM=df_AM_MDM,
      AW=df_AW,
      sherman=df_sherman,
      galagan=df_galagan,
      rustad=df_rustad,
      GSE167232=df_GSE167232,
      )
  return sample_data

def _getGalaganData(is_trinary=True, is_display_errors=False):
  """
  Constructs trinary values for Galagan data.
  These data are normalized and in log2 units.
  The 10h data are Normoxia.
  :param bool is_trinary: convert to Trinary
  :param bool is_display_errors:
  :return pd.DataFrame:
      columns: genes
      index: time instances
      values: trinary based on log2
  """
  df_galagan = transform_data.readGeneCSV(
      csv_file=FILE_GALAGAN)
  dfs = []
  for idx in range(1, 4):
    stg = "rep%d" % idx
    columns = [c for c in df_galagan.columns
        if stg in c]
    col_ref = columns[0]
    columns.remove(col_ref)
    df = df_galagan[columns].copy()
    df = df.apply(lambda c: c - df_galagan[col_ref])
    dfs.append(df)
  df_merge = pd.concat(dfs, axis=1)
  if is_trinary:
    df_trinary = df_merge.applymap(lambda v:
      1 if v >= 1 else -1 if v <= -1 else 0)
    return df_trinary.T
  else:
    return df_merge.T

def serializeFeatureMatrix(df_X, path):
  """
  Serializes the feature vector as a CSV file.
  :param pd.DataFrame df_X:
  :param str path: output file path
  """
  df = df_X.copy()
  df.index.name = cn.INDEX
  df.to_csv(path)

def mkFeatureMatrices(sample_data,
    directory=cn.TRINARY_SAMPLES_DIR):
  """
  Creates data in trinary feature matrix.
  :param SampleData sample_data:
  """
  for source in SAMPLES:
    path = os.path.join(directory, "%s.csv" % source)
    serializeFeatureMatrix(
        sample_data.__getattribute__(source), path)
    

################## CLASSES ###############
class NormalizedData(object):
  """ Exposes values described above. """

  def __init__(self, is_display_errors=True,
      is_averaged=True, is_regulator=False):
    """
    :param bool is_display_errors: Shows errors encountered
    :param bool is_averaged: Use averaged read counts
    :param bool is_regulator: use regulators for TRN
    Public instance variables:
      df_X are normalized read counts
      ser_y - numeric value of state corresponding to each row in df_X
      self.state_dct:
          key: state name
          value: state index
      self.features: list of names of gene
    """
    self._is_display_errors = is_display_errors
    self.provider = DataProvider(
        is_display_errors=self._is_display_errors)
    self.provider.do()
    if is_averaged:
      self.df_X = self.provider.df_normalized.T
    else:
      # Use the adjusted values for each replication
      dfs = [df.copy() for df in
          self.provider.dfs_adjusted_read_count_wrtT0_log2]
      self.df_X = pd.concat([df.T for df in dfs])
    drop_indices = self._getDropIndices(self.df_X.index)
    self.df_X = self.df_X.drop(drop_indices)
    if is_regulator:
      self.df_X = _subsetToRegulators(self.df_X)
    self.features = self.df_X.columns.tolist()
    # Create class information
    ser_y = self.provider.df_stage_matrix[cn.STAGE_NAME]
    if not is_averaged:
      # Replica information has a specical time format
      num_repl = len(self.provider.dfs_read_count)
      sers = []
      for idx in range(num_repl):
        new_ser_y = ser_y.copy()
        new_ser_y.index = self.provider.makeTimes(
           suffix=idx)
        sers.append(new_ser_y)
      ser_y = pd.concat(sers)
    ser_y = ser_y.drop(self._getDropIndices(ser_y.index))
    # Equate Normoxia and Resuscitation if there are too
    # few states
    states = ["Normoxia", "Transition", "Stage1a", "Stage1b", "StageII",
        "Resuscitation"]
    if len(ser_y[ser_y == cn.STATE_NORMOXIA])  \
        <= MIN_NUM_NORMOXIA:
      ser_y[ser_y == cn.STATE_NORMOXIA] = cn.STATE_RESCUSCITATION
      states.remove("Normoxia")
    # Create converter from state name to numeric index
    self.state_dct = {k: v for v, k in enumerate(states)}
    self.ser_y = ser_y.apply( lambda k: self.state_dct[k])

  def _getDropIndices(self, indices,
      drop_index=cn.TIME_0):
    """
    Handles dropping time index when replicas are present
    """
    result = []
    for idx in indices:
      splits = idx.split(data_provider.SEPARATOR)
      if splits[0] == drop_index:
        result.append(idx)
    return result
  

class TrinaryData(NormalizedData):

  def __init__(self, is_dropT1=True,
      is_display_errors=True, **kwargs):
    """
    self.df_X
        columns: gene index
        index: times
        values: trinary
    self.ser_y
        index: times
        values: state index
    self.features: list of names of gene groups
    """
    super().__init__(is_display_errors=is_display_errors,
        **kwargs)
    self.df_X = transform_data.aggregateGenes(
        df=self.df_X)
    drop_indices = self._getDropIndices(self.df_X.index)
    self.df_X = self.df_X.drop(drop_indices)
    if is_dropT1:
      t1_indices = self._getDropIndices(self.df_X.index,
          drop_index=T1_INDEX)
      self.df_X = self.df_X.drop(t1_indices)
      self.ser_y = self.ser_y.drop(t1_indices)
    sorted_index = sorted(self.ser_y.index,
        key=lambda v: float(v[1:]))
    self.ser_y = self.ser_y[sorted_index]
    self.features = self.df_X.columns.tolist()
    #self.df_X.columns = range(len(self.features))
    self.df_fstat = None

  def plotFeatureSignificanceByState(self,
      max_sl=0.25, max_rank=50, is_plot=True,
      figsize=(8, 6)):
    """
    Constructs a heatmap of F-statistic significance 
    levels by state.
    :param float threshold_sl: maximum significance level
    :param int max_rank: number of genes on y-axis
        Genes orderd by significance level across states
    """
    if self.df_fstat is None:
      self.df_fstat = util_classifier.makeFstatDF(
          self.df_X, self.ser_y)
    threshold = -np.log10(max_sl)  # Number of zeroes
    df_plot = self.df_fstat.applymap(
        lambda v: np.nan if v < threshold else v)
    df_plot = df_plot.loc[df_plot.index[
        0:(max_rank-1)], :]
    plt.figure(figsize=figsize)
    ax = plt.gca()
    ax.set_xticks(np.arange(len(df_plot.columns))+0.5)
    ax.set_xticklabels(df_plot.columns)
    ax.set_yticks(np.arange(len(df_plot.index))+0.5)
    ax.set_yticklabels(df_plot.index, rotation=0)
    ax.set_xlabel("State")
    heatmap = plt.pcolor(df_plot)
    _ = plt.colorbar(heatmap)
    if is_plot:
      plt.show()
