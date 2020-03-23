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
import common.transform_data as transform_data

import numpy as np
import pandas as pd

T1_INDEX = "T1"
MIN_NUM_NORMOXIA = 2  # Minimum number of normoxia states


class NormalizedData(object):
  """ Exposes values described above. """

  def __init__(self, is_display_errors=True, is_averaged=True):
    """
    :param bool is_display_errors: Shows errors encountered
    :param bool is_averaged: Use averaged read counts
    Public instance variables:
      df_X are normalized read counts
      states_dict - mapping of literal to numeric values of state
      ser_y - numeric value of state corresponding to each row in df_X
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
    drop_indices = self._getTimeIndices(self.df_X.index)
    self.df_X = self.df_X.drop(drop_indices)
    self.features = self.df_X.columns.tolist()
    self.df_X.columns = range(len(self.features))
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
    ser_y = ser_y.drop(self._getTimeIndices(ser_y.index))
    # Equate Normoxia and Resuscitation if there are too
    # few states
    if len(ser_y[ser_y == cn.STATE_NORMOXIA]) <= MIN_NUM_NORMOXIA:
      ser_y[ser_y == cn.STATE_NORMOXIA]  \
          = cn.STATE_RESCUSCITATION
    # Create converter from state name to numeric index
    states = ser_y.unique()
    self.state_dict = {k: v for v, k in enumerate(states)}
    self.ser_y = ser_y.apply(lambda k: self.state_dict[k] )

  def _getTimeIndices(self, indices, time_index=cn.TIME_0):
    try:
      result = [i for i in indices if cn.TIME_0 in i]
    except:
      import pdb; pdb.set_trace()
    return result
  

class TrinaryData(NormalizedData):

  def __init__(self, is_dropT1=True,
      is_display_errors=True, **kwargs):
    """
    self.df_X are trinary values
    """
    super().__init__(is_display_errors=is_display_errors,
        **kwargs)
    self.df_X = transform_data.aggregateGenes(
        df=self.df_X)
    drop_indices = self._getTimeIndices(self.df_X.index)
    self.df_X = self.df_X.drop(drop_indices)
    if is_dropT1:
      t1_indices = self._getTimeIndices(self.df_X.index,
          time_index=T1_INDEX)
      self.df_X = self.df_X.drop(t1_indices)
      self.ser_y = self.ser_y.drop(t1_indices)
    self.df_X.index = sorted(self.df_X.index,
        key=lambda v: float(v[1:]))
    self.features = self.df_X.columns.tolist()
    self.df_X.columns = range(len(self.features))
