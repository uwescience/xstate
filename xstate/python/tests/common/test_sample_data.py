import common.constants as cn
from common import sample_data
from common_python.util.persister import Persister

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = True
IS_PLOT = True

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PERSISTER_PATH = os.path.join(TEST_DIR, "test_sample_data_persister.pcl")
PERSISTER = Persister(PERSISTER_PATH)
if PERSISTER.isExist():
  SAMPLE_DATA = PERSISTER.get()
else:
  SAMPLE_DATA = sample_data.getSampleData()
  SAMPLE_DATA.initialize()
  PERSISTER.set(SAMPLE_DATA)


class TestSampleData(unittest.TestCase):

  def _remove(self):
    for sample_name in sample_data.SAMPLES:
      path = os.path.join(TEST_DIR, "%s.csv" % sample_name)
      if os.path.isfile(path):
        os.remove(path)

  def isValidDataframe(self, df, columns):
    self.assertTrue(isinstance(df, pd.DataFrame))
    diff = set(df.columns).symmetric_difference(columns)
    self.assertEqual(len(diff), 0)
    self.assertGreater(len(df), 0)

  def setUp(self):
    self._remove()
    self.data = sample_data.SampleData()

  def tearDown(self):
    self._remove()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertIsNone(self.data.df_galagan)

  def testInitialize(self):
    if IGNORE_TEST:
      return
    self.data.initialize()
    for sample_name in sample_data.SAMPLES:
      df = self.data.getDataFrame(sample_name)
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
    if IGNORE_TEST:
      return
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

  def testGetTrinaryFromGeneLists(self):
    if IGNORE_TEST:
      return
    df = sample_data.getTrinaryFromGeneLists(
        induced_path=None, repressed_path=None)
    ser = df[df.columns.tolist()[0]]
    self.assertEqual(ser.apply(lambda v: np.abs(v)).sum(),
        0)
    #
    df = sample_data.getTrinaryFromGeneLists()
    ser = df[df.columns.tolist()[0]]
    self.assertGreater(len(ser[ser == -1]), 0)
    self.assertGreater(len(ser[ser == 1]), 0)

  def testSerialize(self):
    if IGNORE_TEST:
      return
    SAMPLE_DATA.serialize(directory=TEST_DIR)
    for sample_name in sample_data.SAMPLES:
      path = os.path.join(TEST_DIR, "%s.csv" % sample_name)
      self.assertTrue(os.path.isfile(path))

  def testAverageReplicas(self):
    # TESTING
    replica_names = ["d1", "d2", "d3", "d5", "d7", "d8"]
    df = SAMPLE_DATA.df_galagan
    df_result = sample_data.SampleData.averageReplicas(
        SAMPLE_DATA.df_galagan, replica_names)
    diff = set(df.columns).symmetric_difference(df_result.columns)
    self.assertEqual(len(diff), 0)
    diff = set(df_result.index).symmetric_difference(replica_names)
    self.assertEqual(len(diff), 0)
    #
    COL_A = "a"
    COL_B = "b"
    SIZE = 3
    df = pd.DataFrame({
        COL_A: list(range(2*SIZE)),
        COL_B: [10*v for v in range(2*SIZE)]
    })
    indices = ["c1_%d" % n for n in range(SIZE)]
    other_indices = [i.replace("c1", "c2") for i in indices]
    indices.extend(other_indices)
    df.index = indices
    replica_names = ["c1", "c2"]
    df_result = sample_data.SampleData.averageReplicas(df, replica_names)
    expected = (SIZE-1)*(SIZE-2)/2.0
    self.assertEqual(df_result.loc[replica_names[0], COL_A], expected)
    self.assertEqual(df_result.loc[replica_names[0], COL_B], 10*expected)
    

if __name__ == '__main__':
  unittest.main()
