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
from common import util
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
PROVIDER = DataProvider(is_display_errors=False)
PROVIDER.do()


################## FUNCTIONS ###############
def subsetToRegulators(df):
  regulators = PROVIDER.df_trn_unsigned[cn.TF]
  regulators = list(set(regulators))
  regulator_cols = list(set(df.columns).intersection(regulators))
  for column in df.columns:
    if not column in regulator_cols:
      del df[column]

def convertToTrinary(df, threshold_low=-1, threshold_high=1):
  """
  Converts the dataframe to trinary values using the indicated thresholds.

  Parameters
  ----------
  threshold_low: float
      threshold for low value; set to -1
  threshold_high: float
      threshold for high value; set to +1
  
  Returns
  -------
  pd.DataFrame
  """
  df_trinary = df.applymap(lambda v:
    1 if v >= threshold_high else -1 if v <= threshold_low else 0)
  return df_trinary

def serializeFeatureMatrix(df_X, path):
  """
  Serializes the feature vector as a CSV file.
  :param pd.DataFrame df_X:
  :param str path: output file path
  """
  df = df_X.copy()
  df.index.name = cn.INDEX
  df.to_csv(path)
    

################## CLASSES ###############
class NormalizedData(object):
  """ Exposes values described above. """

  def __init__(self, is_averaged=True, is_regulator=False, **kwargs):
    """
    :param bool is_averaged: Use averaged read counts
    :param bool is_regulator: use regulators for TRN
    :param dict kwargs: options passed to DataProvider

    Public instance variables:
      df_X are normalized read counts
          instances are either times (begin with T) for stage (S)
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
      subsetToRegulators(self.df_X)
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

  def __init__(self, is_dropT1=True, is_stage_averaged=False, **kwargs):
    """
    self.df_X
        columns: gene index
        index: times
        values: trinary
    self.ser_y
        index: times
        values: state index
    self.features: list of names of gene groups
    
    
    Parameters
    ----------
    is_dropT1: bool
        Drop time T1
    is_stage_averaged: bool
        Average log units by stage
    
    Returns
    -------
    """
    super().__init__(**kwargs)
    if is_stage_averaged:
      COL_Y = "col_y"
      df_X = self.df_X.copy()
      df_X[COL_Y] = self.ser_y
      self.df_X = df_X.groupby(COL_Y).mean()
      ser_y = pd.Series([n for n in range(len(set(self.ser_y.values)))])
      self.ser_y = ser_y
      indices = ["S%d" % n for n in self.ser_y]
      self.df_X.index = indices
      self.ser_y.index = indices
      # Variability of the aggregation by stage
      self.df_std = df_X.groupby(COL_Y).std()
    else:
      self.df_X = transform_data.aggregateGenes(
          df=self.df_X)
      self.df_std = None
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
      figsize=(8, 6), is_color_bar=True):
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
    if is_color_bar:
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

  def plotExpressionLevels(self, features, df_X=None, ser_y=None,
      is_plot=True, title="", figsize=(20, 5), is_color_bar=True):
    """
    Heat map of expression levels for features. Shades states.
    
    Parameters
    ----------
    features: list-str
    df_X: DataFrame (feature vector)
        if non-None, then this feature vector is used instead of self.df_X
    ser_y: Series (classes)
    """
    if df_X is None:
      df_X = self.df_X
      ser_y = self.ser_y
    # Internal constants
    ROTATION = 30
    FONTSIZE = 14
    # Shade replications
    fig, ax = plt.subplots(1, figsize=figsize)
    columns = list(set(features).intersection(df_X.columns))
    columns.sort()
    new_df_X = df_X.copy()
    new_df_X = new_df_X[columns]
    sns.heatmap(new_df_X.T, cmap="seismic", ax=ax, vmin=-1, vmax=1,
        cbar=is_color_bar)
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
