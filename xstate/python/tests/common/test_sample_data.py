import common.constants as cn
from common import sample_data

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = True
IS_PLOT = True

  


class TestSampleData(unittest.TestCase):

  def isValidDataframe(self, df, columns):
    self.assertTrue(isinstance(df, pd.DataFrame))
    diff = set(df.columns).symmetric_difference(columns)
    self.assertEqual(len(diff), 0)
    self.assertGreater(len(df), 0)

  def setUp(self):
    self.data = sample_data.SampleData()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertIsNone(self.data.df_galagan)

  def testInitialize(self):
    if IGNORE_TEST:
      return
    self.data.initialize()
    for attr in sample_data.SAMPLES:
      df = self.data.__getattribute__(attr)
      self.assertTrue(isinstance(df, pd.DataFrame))

  def testGetGalaganData(self):
    if IGNORE_TEST:
      return
    df = self.data._getGalaganData()
    trues = ["Rv" in c for c in df.columns]
    self.isValidDataframe(df, df.columns)

  def testGetGSE167232(self):
    if IGNORE_TEST:
      return
    is_regulator = False
    is_display_errors = True
    lastDF = self.data._getGSE167232()
    for _ in range(5):
      newDF = self.data._getGSE167232()
      self.assertTrue(newDF.equals(lastDF))

  def testGetSampleData(self):
    # TESTING
    def getDFS(sample):
      return [sample.df_AW, sample.df_AM_MDM, sample.df_galagan,
          sample.df_sherman, sample.df_GSE167232]
    #
    def test_single(**kwargs):
      sample = sample_data.getSampleData(**kwargs)
      for df in getDFS(sample):
        self.isValidDataframe(df, df.columns)
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
        sample_data.REF_TYPE_BIOREACTOR,
        sample_data.REF_TYPE_SELF,
        sample_data.REF_TYPE_POOLED,
        ]
    for ref_type in ref_types:
      sample_large = test_single(is_regulator=False, ref_type=ref_type)
      sample_small = test_single(is_regulator=True, ref_type=ref_type)
      test_greater(sample_large, sample_small)
    

if __name__ == '__main__':
  unittest.main()
