"""Normalized and Trinary Data."""

"""
df_X - feature matrix values with numeric column names.
       T0 is deleted
ser_y - state values. States are numbered to be consistent with the
        order in cn.STATE_NAMES
features - names of genes
state_dct - dictionary mapping string names to numbers

Note that the terms "stage" and "state" are used interchangably.
"""

import common.constants as cn
from common.msg import writeMessage
from common.data_provider import DataProvider
from common import data_provider
import common.transform_data as transform_data
from common_python.classifier import util_classifier

import collections
import copy
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

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
PROVIDER = DataProvider(is_display_errors=False)
PROVIDER.do()
# Reference types
REF_TYPE_BIOREACTOR = "ref_type_bioreactor"  # Use bioreactor data as reference
REF_TYPE_SELF = "ref_type_self"  # Internal reference
REF_TYPE_POOLED = "ref_type_pooled"  # Pool the data to get a reference


################## FUNCTIONS ###############
def _subsetToRegulators(df):
  regulators = PROVIDER.df_trn_unsigned[cn.TF]
  regulators = list(set(regulators))
  regulator_cols = list(set(df.columns).intersection(regulators))
  for column in df.columns:
    if not column in regulator_cols:
      del df[column]
#
def _getTrinaryFromGeneLists(
    repressed_path=SHERMAN_REPRESSED_PATH,
    induced_path=SHERMAN_INDUCED_PATH):
  """
  Creates a feature vector from a list of induced
  and repressed genes.

  Parameters
  ----------
  repressed_path: str
  induced_path: str

  Return
  ------
  pd.DataFrame with single column
      index: gene
      value: trinary     
  """
  genes = PROVIDER.dfs_read_count[0].index.tolist()
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
    is_display_errors=False,
    ref_type=REF_TYPE_BIOREACTOR,
    ):
  """
  Acquires data obtain from other soruces.
  :param bool is_regulator: use regulators for TRN
  :param str ref_type: indicates the type of reference value
      T0 data as the reference data
  :return SampleData. Each element is pd.DataFrame:
      columns: Gene
      index: condition
      values: trinary
  """
  COL_TIME = "time"
  COL_REP = "rep"
  COL_CONDITION = "condition"
  def makeSamplesWithPooledReference(csv_file, calcRef, is_time_columns,
      df_data=None):
    """
    Constructs samples with respect to reference instances within the same
    data set.

    Parameters
    ----------
    csv_file: str
        CSV file for the sample data if df_data == None
    calcRef: Function
        parameters: DataFrame
        returns: Series
    is_time_columns: bool
        has a "time" column
    df_data: dataframe
        count (not log2) to be converted into TrinaryData
    Returns
    -------
    DataFrame: Trinary values of samples
    """
    if df_data is None:
      df_data = transform_data.readGeneCSV(csv_dir=cn.SAMPLES_DIR,
          csv_file=csv_file).T
    df_normalized = PROVIDER.normalizeReadsDF(
        df_data.T, is_time_columns=False).T
    df_normalized = df_normalized.applymap(lambda v: max(v, cn.MIN_VALUE))
    ser_ref = calcRef(df_normalized)
    ser_ref = ser_ref.map(lambda v: max(v, cn.MIN_VALUE))  # No 0's
    ser_ref_log = transform_data.convertToLog2(ser_ref)
    df_normalized_log = transform_data.convertToLog2(df_normalized)
    df_sample_trinary = transform_data.calcTrinaryComparison(
        df_normalized_log.T,
        ser_ref=ser_ref_log, is_convert_log2=False).T
    return df_sample_trinary
  #
  def calcRefFromIndices(df, sel_ref_func):
    ref_idxs = [i for i in df.index if ref_sel_func(i)]
    df_ref = df.loc[ref_idxs, :]
    return df_ref.mean()
  #
  def calcRefPooled(df):
    return df.mean(axis=0)
  #
  # AM/MDM
  if (ref_type == REF_TYPE_BIOREACTOR):
    df_AM_MDM = transform_data.trinaryReadsDF(
        is_display_errors=is_display_errors,
        csv_file=FILE_AM_MDM,
        is_time_columns=False).T
  elif (ref_type == REF_TYPE_SELF):
    ref_sel_func = lambda i: ("AM" in i) and (not "1" in i)
    def calcRef(df):
      return calcRefFromIndices(df, ref_sel_func)
    df_AM_MDM = makeSamplesWithPooledReference(FILE_AM_MDM, calcRef, False)
  else:
    df_AM_MDM = makeSamplesWithPooledReference(FILE_AM_MDM, calcRefPooled, False)
  # AW
  if (ref_type == REF_TYPE_BIOREACTOR):
    df_AW = transform_data.trinaryReadsDF(
        csv_file=FILE_AW,
        is_display_errors=is_display_errors,
        is_time_columns=False).T
  elif (ref_type == REF_TYPE_SELF):
    ref_sel_func = lambda i: ("neg" in i) and (not "1" in i)
    def calcRef(df):
      return calcRefFromIndices(df, ref_sel_func)
    df_AW = makeSamplesWithPooledReference(FILE_AW, calcRef, False)
  else: # Pooled
    df_AW = makeSamplesWithPooledReference(FILE_AW, calcRefPooled, False)
  df_AW = df_AW.sort_index()
  # Galagn data
  if (ref_type == REF_TYPE_BIOREACTOR):
    df_galagan = _getGalaganData(
        is_display_errors=is_display_errors, is_trinary=True)
  elif (ref_type == REF_TYPE_SELF):
    df_data = _getGalaganData(
        is_display_errors=is_display_errors, is_trinary=False)
    ref_sel_func = lambda i: ("d1." in i) and ("rep1" not in i)
    def calcRef(df):
      return calcRefFromIndices(df, ref_sel_func)
    df_galagan = makeSamplesWithPooledReference(None, calcRef, False, df_data=df_data)
  else:  # Pooled
    df_data = _getGalaganData(
        is_display_errors=is_display_errors, is_trinary=False)
    df_galagan = makeSamplesWithPooledReference(None, calcRefPooled, False, df_data=df_data)
  #
  df_sherman = _getTrinaryFromGeneLists()
  df_sherman = df_sherman.transpose()
  # Rustad
  if (ref_type == REF_TYPE_BIOREACTOR):
    df_rustad = transform_data.trinaryReadsDF(
        csv_file=FILE_RUSTAD,
        is_display_errors=is_display_errors,
        is_convert_log2=False,
        is_time_columns=False).T
  elif (ref_type == REF_TYPE_SELF):
    ref_sel_func = lambda i: ("_4hr_" in i) and ("rep6" not in i)
    def calcRef(df):
      return calcRefFromIndices(df, ref_sel_func)
    df_rustad = makeSamplesWithPooledReference(FILE_RUSTAD, calcRef, False)
  else: # pooled
    df_rustad = makeSamplesWithPooledReference(FILE_RUSTAD, calcRefPooled, False)
  # Construct the major sor index for rustad
  time_vals = ["4hr", "8hr", "12hr", "1day", "4day", "7day"]
  reps  = [i.split("_")[-1] for i in df_rustad.index]
  times  = [i.split("_")[-2] for i in df_rustad.index]
  conditions = ["_".join(i.split("_")[0:2]) for i in df_rustad.index]
  df_rustad[COL_TIME] = [time_vals.index(v) for v in times]
  df_rustad[COL_REP] = reps
  df_rustad[COL_CONDITION] = conditions
  df_rustad = df_rustad.sort_values([COL_CONDITION, COL_TIME, COL_REP])
  for col in [COL_REP, COL_TIME, COL_CONDITION]:
    del df_rustad[col]
  # GSE167232
  if (ref_type == REF_TYPE_BIOREACTOR) or (ref_type == REF_TYPE_SELF):
    if (ref_type == REF_TYPE_SELF):
      message = "\n**No self reference defined for GSE167232."
      message += " Using bioreactor data.\n"
      writeMessage(message)
    df_GSE167232 = transform_data.trinaryReadsDF(
        csv_file=FILE_GSE167232,
        is_display_errors=is_display_errors,
        is_normalized=True,
        is_time_columns=False).T
  else: # pooled
    df_GSE167232 = makeSamplesWithPooledReference(FILE_GSE167232,
        calcRefPooled, False)
  # Restrict to regulators?
  if is_regulator:
    for df in [df_AM_MDM, df_AW, df_sherman, df_galagan,
        df_rustad, df_GSE167232]:
      _subsetToRegulators(df)
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

def _getGSE167232(is_regulator, is_display_errors):
  df = transform_data.trinaryReadsDF(
      csv_file=FILE_GSE167232,
      is_display_errors=is_display_errors,
      is_normalized=True,
      is_time_columns=False).T
  if is_regulator:
    _subsetToRegulators(df)
  return df

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

  def __init__(self,
      is_averaged=True, is_regulator=False,
      **kwargs):
    """
    :param bool is_averaged: Use averaged read counts
    :param bool is_regulator: use regulators for TRN
    :param dict kwargs: options passed to DataProvider

    Public instance variables:
      df_X are normalized read counts
      ser_y - numeric value of state corresponding to each row in df_X
      self.state_dct:
          key: state name
          value: state index
      self.features: list of names of gene
    """
    self.provider = DataProvider(**kwargs)
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
      _subsetToRegulators(self.df_X)
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
    states = list(cn.STATE_NAMES)
    ser_y = ser_y.drop(self._getDropIndices(ser_y.index))
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

  def __init__(self, is_dropT1=True, **kwargs):
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
    super().__init__(**kwargs)
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

  def subsetToStates(self, states, genes=None):
    """
    Transforms the trinary data as follows:
       df_X - subsets the genes
       ser_y - creates binary classes that selects a subset of states
    
    Parameters
    ----------
    states: list-str
        str: timepoint
    genes: list-str
    
    Returns
    -------
    TrinaryData
    """
    if genes is None:
      genes = list(self.df_X.columns)
    #
    data = copy.deepcopy(self)
    data.df_X = pd.DataFrame(data.df_X[genes])
    if states is not None:
      ser_y = data.ser_y.copy()
      numeric_states = []
      value_dct = {i: 1 if (s in states) else 0 for  i, s
          in enumerate(self.state_dct.keys())}
      new_values = [value_dct[v] for v in data.ser_y]
      data.ser_y = pd.Series(new_values)
      data.ser_y.index = ser_y.index
    return data 

  def getStateNames(self, state_ints):
    """
    Provides the name of the state for the value of a state.
    
    Parameters
    ----------
    states: list-int

    Return
    list-str
    """
    rev_dct = {v: k for k, v in self.state_dct.items()}
    return [rev_dct[i] for i in state_ints]

  def plotExpressionLevels(self, features, df_X=None, is_plot=True, title=""):
    """
    Heat map of expression levels for features. Shades states.
    
    Parameters
    ----------
    features: list-str
    df_X: DataFrame (feature vector)
        if non-None, then this feature vector is used instead of self.df_X
    """
    if df_X is None:
      df_X = self.df_X
      ser_y = self.ser_y
    else:
      ser_y = None
    # Internal constants
    ROTATION = 30
    FONTSIZE = 14
    # Shade replications
    fig, ax = plt.subplots(1, figsize=(20, 5))
    columns = list(set(features).intersection(df_X.columns))
    columns.sort()
    new_df_X = df_X[columns]
    sns.heatmap(new_df_X.T, cmap="seismic", ax=ax, vmin=-1, vmax=1)
    # Shade the classes
    if ser_y is not None:
      alphas = [0.0, 0.4]
      alpha_idx = 0
      indices = list(ser_y.index)
      for idx, val in enumerate(ser_y.values):
        timepoint = indices[idx]
        stage = PROVIDER.getStateNameForTimepoint(timepoint)
        if (idx == 0):
          ax.text(idx, 0, "%s" % stage, fontsize=FONTSIZE, rotation=ROTATION)
        elif ser_y.values[idx-1] != val:
          ax.text(idx, 0, "%s" % stage, fontsize=FONTSIZE, rotation=ROTATION)
          alpha_idx = 1 - alpha_idx
        ax.axvspan(idx, idx+1, facecolor='grey', alpha=alphas[alpha_idx])
    # Other plot characteristics
    ax.set_title(title, fontsize=18, pad=18, loc="right")
    #
    if is_plot:
      plt.show()
    else:
      plt.close()
