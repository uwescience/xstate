import common.constants as cn
from common import trinary_data
from common.trinary_data import TrinaryData, NormalizedData
from common import trinary_data
from common_python.testing import helpers

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = True
IS_PLOT = True
NUM_REPL = 3
DIR = os.path.dirname(os.path.abspath(__file__))
TEST_SAMPLE_PATH = os.path.join(DIR, "sample.csv")


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

  def _remove(self):
    for path in [TEST_SAMPLE_PATH]:
      if os.path.isfile(path):
        os.remove(path)

  def setUp(self):
    self._remove()

  def tearDown(self):
    self._remove()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    for cls in [TrinaryData, NormalizedData]:
      data = cls()
      self.assertTrue(isinstance(data.df_X, pd.DataFrame))
      self.assertTrue(isinstance(data.ser_y, pd.Series))
      self.assertTrue(isinstance(data.features, list))
      self.assertTrue(isinstance(data.state_dct, dict))
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

  def testGetGalaganData(self):
    if IGNORE_TEST:
      return
    df = trinary_data._getGalaganData(False)
    trues = ["Rv" in c for c in df.columns]
    self.assertTrue(helpers.isValidDataFrame(df,
        df.columns))

  def testGetSampleData(self):
    # TESTING
    def getDFS(sample):
      return [sample.AW, sample.AM_MDM, sample.galagan,
          sample.sherman]
    #
    def test_single(is_regulator=False, is_bioreactor_ref=True):
      sample = trinary_data.getSampleData(
          is_display_errors=False,
          is_bioreactor_ref=is_bioreactor_ref,
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
    # Non-default reference
    sample_not_bioreactor_ref = test_single(is_regulator=False,
        is_bioreactor_ref=False)
    sample_reg = test_single(is_regulator=True)
    sample_full = test_single(is_regulator=False)
    test_greater(sample_full, sample_reg)
    for name in trinary_data.SAMPLES:
      df_not_bioreactor = sample_not_bioreactor_ref.__getattribute__(name)
      df_bioreactor = sample_reg.__getattribute__(name)
      if name not in ["sherman", "GSE167232"]:
        self.assertGreater(len(df_bioreactor), len(df_not_bioreactor))
      else:
        self.assertEqual(len(df_bioreactor), len(df_not_bioreactor))

  def testGetTrinaryFromGeneLists(self):
    if IGNORE_TEST:
      return
    df = trinary_data._getTrinaryFromGeneLists(
        induced_path=None, repressed_path=None)
    ser = df[df.columns.tolist()[0]]
    self.assertEqual(ser.apply(lambda v: np.abs(v)).sum(),
        0)
    #
    df = trinary_data._getTrinaryFromGeneLists()
    ser = df[df.columns.tolist()[0]]
    self.assertGreater(len(ser[ser == -1]), 0)
    self.assertGreater(len(ser[ser == 1]), 0)

  def testSerializeFeatureMatrix(self):
    if IGNORE_TEST:
      return
    sample_data = trinary_data.getSampleData()
    df_X = sample_data.galagan
    self.assertFalse(os.path.isfile(TEST_SAMPLE_PATH))
    trinary_data.serializeFeatureMatrix(df_X,
        TEST_SAMPLE_PATH)
    self.assertTrue(os.path.isfile(TEST_SAMPLE_PATH))

  def testSerializeFeatureMatrix(self):
    if IGNORE_TEST:
      return
    def removeFiles():
      for source in trinary_data.SAMPLES:
        path = os.path.join(DIR, "%s.csv" % source)
        if os.path.isfile(path):
          os.remove(path)
    #
    removeFiles()
    sample_data = trinary_data.getSampleData()
    trinary_data.mkFeatureMatrices(sample_data,
        directory=DIR)
    for source in trinary_data.SAMPLES:
      path = os.path.join(DIR, "%s.csv" % source)
      self.assertTrue(os.path.isfile(path))
    removeFiles()
 
    

if __name__ == '__main__':
  unittest.main()
