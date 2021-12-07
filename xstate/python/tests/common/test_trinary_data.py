import common.constants as cn
from common import trinary_data
from common.data_provider import DataProvider
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
PROVIDER = DataProvider()
PROVIDER.do()


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

  def testTrinaryRefPooled(self):
    if IGNORE_TEST:
      return
    trinary1 = TrinaryData()
    trinary2 = TrinaryData(calcRef=PROVIDER.calcRefPooled)
    self.assertFalse(trinary1.df_X.equals(trinary2.df_X))

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

  def testGetGSE167232(self):
    if IGNORE_TEST:
      return
    is_regulator = False
    is_display_errors = True
    lastDF = trinary_data._getGSE167232(is_regulator, is_display_errors)
    for _ in range(5):
      newDF = trinary_data._getGSE167232(is_regulator, is_display_errors)
      self.assertTrue(newDF.equals(lastDF))

  def testGetSampleData(self):
    if IGNORE_TEST:
      return
    def getDFS(sample):
      return [sample.AW, sample.AM_MDM, sample.galagan,
          sample.sherman, sample.GSE167232]
    #
    def test_single(**kwargs):
      sample = trinary_data.getSampleData(**kwargs)
      for df in getDFS(sample):
        self.assertTrue(helpers.isValidDataFrame(df, df.columns))
        self.assertTrue(np.abs(df.mean().mean()) < np.inf) 
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
    ref_types = [
        trinary_data.REF_TYPE_BIOREACTOR,
        trinary_data.REF_TYPE_SELF,
        trinary_data.REF_TYPE_POOLED,
        ]
    for ref_type in ref_types:
      sample_large = test_single(is_regulator=False, ref_type=ref_type)
      sample_small = test_single(is_regulator=True, ref_type=ref_type)
      test_greater(sample_large, sample_small)

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

  def testSubsetToStages(self):
    # TESTING
    trinary = TrinaryData()
    GENES = ["Rv1927", "Rv3083"]
    subset_trinary = trinary.subsetToStages(["Transition"], genes=GENES)
    len1 = len(trinary.ser_y[trinary.ser_y > 0])
    len2 = len(subset_trinary.ser_y[subset_trinary.ser_y > 0])
    self.assertGreater(len1, len2)
    self.assertEqual(len(GENES), len(subset_trinary.df_X.columns))
    

if __name__ == '__main__':
  unittest.main()
