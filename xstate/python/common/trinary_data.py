"""Normalized and Trinary Data."""

"""
df_X - feature matrix values with numeric column names.
       T0 is deleted
ser_y - state values
features - names of genes
state_dict - dictionary mapping string names to numbers
"""

import common.constants as cn
from common.data_provider import DataProvider
from common import data_provider
import common.transform_data as transform_data
from common_python.classifier import util_classifier

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

T1_INDEX = "T1"
MIN_NUM_NORMOXIA = 2  # Minimum number of normoxia states


################## CLASSES ###############
class NormalizedData(object):
  """ Exposes values described above. """

  def __init__(self, is_display_errors=True, is_averaged=True):
    """
    :param bool is_display_errors: Shows errors encountered
    :param bool is_averaged: Use averaged read counts
    Public instance variables:
      df_X are normalized read counts
      ser_y - numeric value of state corresponding to each row in df_X
      self.state_dict:
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
    if len(ser_y[ser_y == cn.STATE_NORMOXIA])  \
        <= MIN_NUM_NORMOXIA:
      ser_y[ser_y == cn.STATE_NORMOXIA]  \
          = cn.STATE_RESCUSCITATION
    # Create converter from state name to numeric index
    states = ser_y.unique()
    self.state_dict = {k: v for v, k in enumerate(states)}
    self.ser_y = ser_y.apply(
        lambda k: self.state_dict[k])

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
