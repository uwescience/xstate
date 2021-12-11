from common import data_provider
from common.trinary_data import TrinaryData
import common.constants as cn
from common_python.testing import helpers
from common_python.util.persister import Persister
import common_python.util.dataframe as dataframe

import copy
import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
SIZE = 4


class TestDataProvider(unittest.TestCase):

  def setUp(self):
    self.provider = data_provider.DataProvider()
 
  def tearDown(self):
    persister = Persister(cn.DATA_PROVIDER_PERSISTER_PATH)
    persister.remove()

  def init(self):
    self.df_data = self.provider._makeDFFromCSV(
        data_provider.FILENAME_READS)
    self.df_gene_description =  \
        self.provider._makeGeneDescriptionDF()
    self.provider.df_gene_description = self.df_gene_description
    self.df_data = self.df_data.set_index(cn.GENE_ID)

  def makeData(self, size=SIZE):
    df_data = pd.DataFrame({'a': range(10)})
    self.provider._dfs_centered_adjusted_read_count = [df_data for _ in range(SIZE)]

  def checkDF(self, df, is_check_index=True,
      is_check_column=True, is_check_times=False,
      **kwargs):
    """
    Verifies DataFrames
    :param pd.DataFrame df:
    :param bool is_check_index: checks that index is GENE
    :param bool is_check_column: checks column for time format
    :param bool is_check_times: checks the time format
    :param dict kwargs: arguments passed to checkTimes
    """
    # Non-zero length
    self.assertGreater(len(df), 0)
    # Has the right index
    if is_check_index:
      b = set(df.index).issubset( self.df_gene_description.index)
      if not b:
        import pdb; pdb.set_trace()
      self.assertTrue(b)
    # No nan values
    types = [np.dtype('int64'), np.dtype('float64'), np.dtype('bool')]
    if is_check_column:
      for column in df.columns:
        ser = df[column]
        if ser.dtype in types:
          is_nan = np.isnan(ser.sum(skipna=False))
          self.assertFalse(is_nan)
    if is_check_times:
      self.checkTimes(df.columns, **kwargs)

  def checkTimes(self, times, is_replicated=False):
    """
    Verifies that times have the correct format.
    :param list times:
    :param bool is_replicated: expects a .%d format
    """
    columns = []
    for time in times:
      self.assertTrue("T" in time)
      if is_replicated:
        self.assertTrue(data_provider.SEPARATOR in time)
        splits = time.split(data_provider.SEPARATOR)
        columns.append(splits[0])
      else:
        columns.append(time)
    diff = set(columns).symmetric_difference(
        self.provider.df_normalized.columns)
    self.assertEqual(len(diff), 0)
  
  def testEquals(self):
    if IGNORE_TEST:
      return
    self.assertTrue(self.provider.equals(self.provider))
    provider = copy.deepcopy(self.provider)
    self.provider.do()
    self.assertFalse(self.provider.equals(provider))
    self.assertTrue(self.provider.equals(self.provider))

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.init()
    self.assertTrue("data" in self.provider._data_dir)

  def testmakeDFFromCSV(self):
    if IGNORE_TEST:
      return
    self.init()
    df = self.provider._makeDFFromCSV(data_provider.FILENAME_HYPOXIA)
    self.assertGreater(len(df), 0)

  def testMakeHypoxiaDF(self):
    if IGNORE_TEST:
      return
    self.init()
    df = self.provider._makeHypoxiaDF()
    expected_columns = [
        cn.STD, cn.CV, cn.HOURS, cn.SAMPLE, cn.MEAN, 0, 1, 2]
    trues = [c in df.columns for c in expected_columns]
    self.assertTrue(all(trues))
    self.assertGreater(len(df), 0)

  def testMakeGeneDescriptionDF(self):
    if IGNORE_TEST:
      return
    self.init()
    df = self.provider._makeGeneDescriptionDF()
    self.checkDF(df, is_check_column=False)

  def testGetNumRepl(self):
    if IGNORE_TEST:
      return
    self.init()
    self.provider._dfs_centered_adjusted_read_count = range(3)
    self.assertEqual(self.provider._getNumRepl(), 3)

  def testMakeMeanDF(self):
    if IGNORE_TEST:
      return
    self.init()
    self.makeData()
    df = self.provider._makeMeanDF()
    df = df.applymap(lambda v: int(v))
    self.assertTrue(df.equals(self.provider.dfs_centered_adjusted_read_count[0]))

  def testMakeStdDF(self):
    if IGNORE_TEST:
      return
    self.init()
    self.makeData()
    df = self.provider._makeStdDF()
    self.assertTrue(np.isclose(df.sum().sum(), 0))

  def testReduceDF(self):
    if IGNORE_TEST:
      return
    self.init()
    df = self.provider._reduceDF(self.df_data)
    self.assertGreater(len(self.df_data), len(df))
    difference = set(df.columns).symmetric_difference(
        self.df_data.columns)
    self.assertEqual(len(difference), 0)

  def concatDFS(self):
    dfs = []
    dfs.extend(self.provider.dfs_read_count)
    dfs.extend(self.provider.dfs_adjusted_read_count)
    dfs.extend(self.provider.dfs_adjusted_read_count_wrtT0)
    dfs.extend(
        self.provider.dfs_adjusted_read_count_wrtT0_log2)
    dfs.extend(
        self.provider.dfs_centered_adjusted_read_count)
    return dfs

  def testDo(self):
    if IGNORE_TEST:
      return
    def testLessEqual(dfs1, dfs2):
      for idx in range(len(dfs1)):
        df1 = dfs1[idx]
        df2 = dfs2[idx]
        for gene in df1.index:
          ser = df1.loc[gene, :] <= df2.loc[gene, :]
          self.assertEqual(ser.sum(),  len(ser))
    #
    self.init()
    self.provider.do()
    # Specific tests
    dfs_adjusted_read_count  \
        = self.provider.dfs_adjusted_read_count
    dfs_adjusted_read_count_wrtT0  \
        = self.provider.dfs_adjusted_read_count_wrtT0
    dfs_adjusted_read_count_wrtT0_log2  \
        = self.provider.dfs_adjusted_read_count_wrtT0_log2
    dfs_centered_adjusted_read_count  \
        = self.provider.dfs_centered_adjusted_read_count
    testLessEqual(dfs_centered_adjusted_read_count,
        dfs_adjusted_read_count)
    # Lists
    self.assertTrue(isinstance(self.provider.tfs, list))
    self.assertGreater(len(self.provider.tfs), 0)
    # Common tests
    self.assertEqual(
        len(dfs_centered_adjusted_read_count),
        data_provider.NUM_REPL)
    [self.checkDF(df, is_replicated=True) for df in 
        dfs_centered_adjusted_read_count]
    dfs = [
        self.provider.df_gene_description,
        self.provider.df_mean,
        self.provider.df_std,
        self.provider.df_normalized,
        self.provider.df_gene_expression_state,
        ]
    for idx, df in enumerate(dfs):
      self.checkDF(df, is_check_index=False,
          is_check_times=False)
    for idx, df in enumerate(self.concatDFS()):
      self.checkDF(df, is_replicated=True,
          is_check_times=True)
    dfs = [
        self.provider.df_cv,
        self.provider.df_kegg_gene_pathways, 
        self.provider.df_go_terms, 
        self.provider.df_ec_terms, 
        self.provider.df_ko_terms, 
        self.provider.df_kegg_pathways,
        self.provider.df_trn_signed,
        self.provider.df_trn_unsigned,
        ]
    for idx, df in enumerate(dfs):
      self.checkDF(df, is_check_index=False,
        is_check_column=False)
    columns = self.provider.df_stage_matrix.columns
    diff = set(columns).symmetric_difference(
        [cn.STAGE_NAME, cn.STAGE_COLOR])
    self.assertEqual(len(diff), 0)
    self.assertGreater(len(self.provider.df_stage_matrix), 0)

  def testPersistence(self):
    if IGNORE_TEST:
      return
    self.provider.do()
    provider = data_provider.DataProvider()
    provider.do()
    self.provider.equals(provider)

  def testNormalizeReadsDF(self):
    if IGNORE_TEST:
      return
    provider = data_provider.DataProvider(
        is_only_qgenes=False, is_display_errors=False)
    provider.do()
    df = provider.dfs_read_count[0]
    df_normalized = provider.normalizeReadsDF(df)
    columns = ["T%d" % n for n in range(len(df.columns))]
    self.assertTrue(helpers.isValidDataFrame(df_normalized,
        columns))
    self.assertEqual(len(df), len(df_normalized))
    #ser_length = provider.df_gene_description[cn.LENGTH]

  def testGetStates(self):
    if IGNORE_TEST:
      return
    self.provider.do()
    def test(timepoints):
      results = self.provider.getStates(timepoints)
      self.assertEqual(results[0], "Normoxia")
      if len(results) > 1:
        self.assertEqual(results[1], "Resuscitation")
    #
    test(["T1", "T25"])
    test(["T1"])
    test([1, 25])
    test([1])

  def testGetStateNames(self):
    if IGNORE_TEST:
      return
    self.provider.do()
    trinary = TrinaryData()
    result1s = self.provider.getStateNames(trinary.ser_y)
    count = len(set(trinary.ser_y.values))
    self.assertEqual(len(result1s), count)
    #
    ser_y = trinary.ser_y.copy()
    indices = [i + ".0" for i in ser_y.index]
    ser_y.index = indices
    result2s = self.provider.getStateNames(trinary.ser_y)
    self.assertTrue(all([v1 == v2 for v1, v2 in zip (result1s, result2s)]))

  def testCalcRefPooled(self):
    if IGNORE_TEST:
      return
    df = self.provider._getLog2NormalizedReadcounts()
    ser = self.provider.calcRefPooled(df)
    trues = [v >= 0 for v in ser.values]
    self.assertTrue(all(trues))

  def testMakeNormalizedDF(self):
    if IGNORE_TEST:
      return
    provider = data_provider.DataProvider()
    self.provider.do()
    df1 = self.provider._makeNormalizedDF()
    provider = data_provider.DataProvider(calcRef=self.provider.calcRefPooled)
    provider.do()
    df2 = provider._makeNormalizedDF()
    self.assertTrue(isinstance(df1, pd.DataFrame))
    self.assertTrue(isinstance(df2, pd.DataFrame))
    diff = set(df1.columns).symmetric_difference(df2.columns)
    self.assertEqual(len(diff), 0)
    self.assertEqual(len(df1), len(df2))
    self.assertFalse(df1.equals(df2))

  def testGetStateNameForTimepoint(self):
    if IGNORE_TEST:
      return
    self.provider.do()
    name = self.provider.getStateNameForTimepoint("T0")
    self.assertEqual(name, "Normoxia")
    name = self.provider.getStateNameForTimepoint("T24")
    self.assertEqual(name, "Resuscitation")
 


if __name__ == '__main__':
  unittest.main(failfast=True)
