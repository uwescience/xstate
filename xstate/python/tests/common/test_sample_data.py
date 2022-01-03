import common.constants as cn
from common import sample_data
from common_python.util.persister import Persister

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PERSISTER_PATH = os.path.join(TEST_DIR, "test_sample_data_persister.pcl")
PERSISTER = Persister(PERSISTER_PATH)
got_sample = False
if PERSISTER.isExist():
  SAMPLE_DATA = PERSISTER.get()
  if SAMPLE_DATA is None:
    got_sample = False
  else:
    got_sample = True
if not got_sample:
  try:
    SAMPLE_DATA = sample_data.getSampleData()
    SAMPLE_DATA.initialize()
  except:
    SAMPLE_DATA = None
    print("***Proceeding without SAMPLE_DATA")
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
    self.assertGreater(len(df.columns), 0)

  def setUp(self):
    self._remove()
    self.data = sample_data.SampleData()

  def tearDown(self):
    self._remove()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertIsNone(self.data.df_galagan)

  def testInitializeBase(self):
    if IGNORE_TEST:
      return
    data_base = sample_data.SampleData(is_regulator=True)
    data_base.initialize()
    for sample_name in sample_data.SAMPLES:
      df = data_base.getDataframe(sample_name)
      self.assertTrue(isinstance(df, pd.DataFrame))
    # Regulators
    data_regulator = sample_data.SampleData(is_regulator=True)
    data_regulator.initialize()
    for sample_name in sample_data.SAMPLES:
      df = data_regulator.getDataframe(sample_name)
      self.assertTrue(isinstance(df, pd.DataFrame))
      df_base = data_base.getDataframe(sample_name)
      self.assertGreaterEqual(len(df_base.columns), len(df.columns))
      self.assertEqual(len(df_base), len(df))
    # Average
    data_average = sample_data.SampleData(is_average=True)
    data_average.initialize()
    for sample_name in sample_data.SAMPLES:
      df = data_average.getDataframe(sample_name)
      self.assertTrue(isinstance(df, pd.DataFrame))
      df_base = data_base.getDataframe(sample_name)
      self.assertEqual(len(df_base.columns), len(df.columns))
      self.assertGreaterEqual(len(df_base), len(df))

  def testInitializeRefTypeSelf(self):
    if IGNORE_TEST:
      return
    data = sample_data.SampleData(ref_type=sample_data.REF_TYPE_SELF)
    data.initialize()
    for sample_name in sample_data.SAMPLES:
      df = data.getDataframe(sample_name)
      if df is not None:
        self.assertTrue(isinstance(df, pd.DataFrame))

  def testInitializeRefTypePooled(self):
    if IGNORE_TEST:
      return
    data = sample_data.SampleData(ref_type=sample_data.REF_TYPE_POOLED)
    data.initialize()
    for sample_name in sample_data.SAMPLES:
      df = data.getDataframe(sample_name)
      self.assertTrue(isinstance(df, pd.DataFrame))

  def testRefTypeSelfRustad(self):
    if IGNORE_TEST:
      return
    sample = sample_data.getSampleData(ref_type=sample_data.REF_TYPE_SELF)
    df = sample["rustad"]
    self.isValidDataframe(df, df.columns)

  def testGetSampleData(self):
    if IGNORE_TEST:
      return
    def getDFS(sample):
      return [sample.df_AW, sample.df_AM_MDM, sample.df_galagan,
          sample.df_GSE167232, sample.df_rustad]
    #
    def test_single(**kwargs):
      sample = sample_data.getSampleData(**kwargs)
      for df in getDFS(sample):
        if df is not None:
          self.isValidDataframe(df, df.columns)
          self.assertTrue(np.abs(df.mean().mean()) < np.inf) 
      return sample
    #
    def test_greater(sample_large, sample_small):
      dfs_large = getDFS(sample_large)
      dfs_small = getDFS(sample_small)
      for idx in range(len(dfs_large)):
        if dfs_large[idx] is None:
          self.assertIsNone(dfs_small[idx])
        else:
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

  def testSerialize(self):
    if IGNORE_TEST:
      return
    SAMPLE_DATA.serialize(directory=TEST_DIR)
    for sample_name in sample_data.SAMPLES:
      path = os.path.join(TEST_DIR, "%s.csv" % sample_name)
      self.assertTrue(os.path.isfile(path))

  def testAverageReplicas(self):
    if IGNORE_TEST:
      return
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

  def testGetItem(self):
    if IGNORE_TEST:
      return
    df_AW = SAMPLE_DATA["AW"]
    self.assertTrue(isinstance(df_AW, pd.DataFrame))

  def testKeys(self):
    if IGNORE_TEST:
      return
    samples = SAMPLE_DATA.keys()
    diff = set(samples).symmetric_difference(sample_data.SAMPLES)
    self.assertEqual(len(diff), 0)

  def testValues(self):
    if IGNORE_TEST:
      return
    dfs = SAMPLE_DATA.values()
    self.assertEqual(len(dfs), len(sample_data.SAMPLES))
    for df in dfs:
      self.assertTrue(isinstance(df, pd.DataFrame))

  def testMakeSortKey(self):
    if IGNORE_TEST:
      return
    result1 = SAMPLE_DATA._makeSortKey("AW", "AW_plus_1")
    result2 = SAMPLE_DATA._makeSortKey("AW", "AW_neg_2")
    self.assertGreater(result2, result1)
    

if __name__ == '__main__':
  unittest.main()
