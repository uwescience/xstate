import common.constants as cn
from common.trinary_data import TrinaryData, NormalizedData
from common import trinary_data
from common_python.testing import helpers

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
NUM_REPL = 3


################### FUNCTIONS ############
def isConsistentState(ser_y):
  # Check consistency of states
  times = ser_y.index
  if len(times) == 0:
    return
  if not isinstance(times[0], str):
    return
  df = pd.DataFrame()
  df["state"] = ser_y
  df["time"] = [t.split(".")[0] for t in ser_y.index]
  df_g = df.groupby("time")
  trues = []
  for time in df_g.groups.keys():
    serr = df[df["time"]==time]["state"]
    if len(serr.unique()) > 1:
      import pdb; pdb.set_trace()
      pass


################### Tests ############
class TestTrinaryData(unittest.TestCase):

  def setUp(self):
    pass

  def testConstructor(self):
    if IGNORE_TEST:
      return
    for cls in [TrinaryData, NormalizedData]:
      data = cls()
      self.assertTrue(isinstance(data.df_X, pd.DataFrame))
      self.assertTrue(isinstance(data.ser_y, pd.Series))
      self.assertTrue(isinstance(data.features, list))
      self.assertTrue(isinstance(data.state_dict, dict))
      self.assertTrue(helpers.isValidDataFrame(data.df_X,
          data.df_X.columns))
      self.assertEqual(len(data.df_X.columns),
          len(data.features))
      self.assertEqual(len(data.df_X),
          len(data.ser_y))

  def testNonAveraged(self):
    if IGNORE_TEST:
      return
    def test(df_X, ser_y):
      isConsistentState(ser_y)
      self.assertEqual(len(df_X), len(ser_y))
    #
    data1 = TrinaryData(is_averaged=False,
        is_dropT1=False)
    data2 = TrinaryData(is_averaged=True,
        is_dropT1=False)
    for data in [data1, data2]:
      test(data.df_X, data.ser_y)
    self.assertGreater(len(data1.df_X), len(data2.df_X))
    self.assertGreater(len(data1.df_X.columns),
        len(data2.df_X.columns))
    # Replicated data should have 3 Normoxia states
    self.assertEqual(data1.ser_y[data1.ser_y==0].count(),
          NUM_REPL)

  def testPlotFeatureSignificanceByState(self):
    if IGNORE_TEST:
      return
    trinary = TrinaryData(is_averaged=False,
        is_dropT1=False)
    trinary.plotFeatureSignificanceByState(
        is_plot=IS_PLOT)

  def testRegulator(self):
    if IGNORE_TEST:
      return
    trinary_full = TrinaryData(is_averaged=False,
        is_dropT1=False, is_regulator=False)
    trinary_regulator = TrinaryData(is_averaged=False,
        is_dropT1=False, is_regulator=True)
    self.assertGreater(len(trinary_full.df_X.columns),
        len(trinary_regulator.df_X.columns))

  def testGetSampleData(self):
    if IGNORE_TEST:
      return
    def getDFS(sample):
      return [sample.AW, sample.AM_MDM, sample.galagan]
    #
    def test_single(is_regulator=False):
      sample = trinary_data.getSampleData(
          is_display_errors=False,
          is_regulator=is_regulator)
      for df in getDFS(sample):
        self.assertTrue(helpers.isValidDataFrame(df,
            df.columns))
      return sample
    #
    def test_greater(sample_large, sample_small):
      dfs_large = getDFS(sample_large)
      dfs_small = getDFS(sample_small)
      for idx in range(len(dfs_large)):
        self.assertGreater(
            len(dfs_large[idx].columns),
            len(dfs_small[idx].columns))
    #
    sample_reg = test_single(is_regulator=True)
    sample_full = test_single(is_regulator=False)
    test_greater(sample_full, sample_reg)
    

if __name__ == '__main__':
  unittest.main()
