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
    self.df_X = self.provider.df_normalized.T
    self.df_X = self.df_X.drop(index="T0")
    self.features = self.df_X.columns.tolist()
    self.df_X.columns = range(len(self.features))
    # Create class information
    ser_y = self.provider.df_stage_matrix[cn.STAGE_NAME]
    ser_y = ser_y.drop(index="T0")
    ser_y = ser_y.copy()
    ser_y[ser_y == 'Normoxia'] = 'Resuscitation'
    # Create converter from state name to numeric index
    states = ser_y.unique()
    self.state_dict = {k: v for v, k in enumerate(states)}
    self.ser_y = ser_y.apply(lambda k: self.state_dict[k] )


class TrinaryData(NormalizedData):

  def __init__(self, is_dropT1=True, is_display_errors=True):
    """
    self.df_X are trinary values
    """
    super().__init__(is_display_errors=is_display_errors)
    self.df_X = transform_data.aggregateGenes(provider=self.provider)
    self.df_X = self.df_X.T
    self.df_X = self.df_X.drop(index="T0")
    if is_dropT1:
      self.df_X = self.df_X.drop(index="T1")
      self.ser_y = self.ser_y.drop(index="T1")
    self.features = self.df_X.columns.tolist()
    self.df_X.columns = range(len(self.features))
