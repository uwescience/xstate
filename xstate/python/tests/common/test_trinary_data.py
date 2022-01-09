import common.constants as cn
from common import trinary_data
from common import data_provider
from common.data_provider import DataProvider
from common import sample_data
from common.trinary_data import TrinaryData, NormalizedData
from common import trinary_data
from common_python.testing import helpers
from common_python.util.persister import Persister

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
NUM_REPL = 3
DIR = os.path.dirname(os.path.abspath(__file__))
TEST_SAMPLE_PATH = os.path.join(DIR, "test_trinary_data_sample.csv")
PERSISTER_PATH = os.path.join(DIR, "test_trinary_data_persister.pcl")
PERSISTER = Persister(PERSISTER_PATH)
GENES = ["Rv1927", "Rv3083"]
if PERSISTER.isExist():
  PROVIDER, SAMPLE_DATA = PERSISTER.get()
else:
  SAMPLE_DATA = sample_data.getSampleData()
  PROVIDER = DataProvider(is_reinitialize=True)
  PROVIDER.do()
  PERSISTER.set((PROVIDER, SAMPLE_DATA))


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
    trinary1 = TrinaryData(is_reinitialize=True)
    trinary2 = TrinaryData(calcRef=PROVIDER.calcRefPooled, is_reinitialize=True)
    is_different = False
    for column in trinary1.df_X.columns:
      if column not in trinary2.df_X.columns:
        is_different = True
        break
      if not trinary1.df_X[column].equals(trinary2.df_X[column]):
        is_different = True
        break
    self.assertTrue(is_different)

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

  def testSubsetToStates(self):
    if IGNORE_TEST:
      return
    trinary = TrinaryData()
    subset_trinary = trinary.subsetToStates(["Transition"], genes=GENES)
    len1 = len(trinary.ser_y[trinary.ser_y > 0])
    len2 = len(subset_trinary.ser_y[subset_trinary.ser_y > 0])
    self.assertGreater(len1, len2)
    self.assertEqual(len(GENES), len(subset_trinary.df_X.columns))

  def testGetStateNames(self):
    if IGNORE_TEST:
      return
    trinary = TrinaryData()
    names = trinary.getStateNames([0, 4])
    self.assertEqual(names[0], "Transition")
    self.assertEqual(names[1], "Resuscitation")

  def testPlotExpressionLevels(self):
    if IGNORE_TEST:
      return
    trinary = TrinaryData()
    trinary.plotExpressionLevels(GENES, is_plot=IS_PLOT, title="title")
    trinary.plotExpressionLevels(GENES, df_X=trinary.df_X, is_plot=IS_PLOT)
    trinary.plotExpressionLevels(GENES, df_X=trinary.df_X, 
        ser_y=trinary.ser_y, is_plot=IS_PLOT)
    trinary.plotExpressionLevels(GENES, is_plot=IS_PLOT, title="title",
        is_color_bar=False)

  def testSerializeFeatureMatrix(self):
    if IGNORE_TEST:
      return
    df_X = SAMPLE_DATA.df_galagan
    self.assertFalse(os.path.isfile(TEST_SAMPLE_PATH))
    trinary_data.serializeFeatureMatrix(df_X,
        TEST_SAMPLE_PATH)
    self.assertTrue(os.path.isfile(TEST_SAMPLE_PATH))
    

if __name__ == '__main__':
  unittest.main()
