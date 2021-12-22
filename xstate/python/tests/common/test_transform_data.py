from common import transform_data
import common.constants as cn
from common.trinary_data import TrinaryData
from common.data_provider import DataProvider
from common_python.testing import helpers

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
SER = pd.Series([float(v) for v in range(10)])
DF = pd.DataFrame({"a": SER})
DF["b"] = DF["a"]*10


class TestFunctions(unittest.TestCase):

  def setUp(self):
    if IGNORE_TEST:
      return
    self._init()

  def _init(self):
    self.provider = DataProvider()
    self.provider.do()

  def testMakeTrinaryData(self):
    if IGNORE_TEST:
      return
    df = transform_data.makeTrinaryData(
        df=self.provider.df_normalized)
    columns = self.provider.df_normalized.columns
    self.assertTrue(helpers.isValidDataFrame(df, columns))

  def testAggregateGenes(self):
    if IGNORE_TEST:
      return
    provider = DataProvider()
    provider.do()
    df = transform_data.aggregateGenes(provider=provider)
    self.assertTrue(helpers.isValidDataFrame(df,
        provider.df_normalized.columns))

  def testTrinaryReadsDF1(self):
    if IGNORE_TEST:
      return
    provider = DataProvider()
    provider.do()
    df = provider.dfs_read_count[0]
    df_result = transform_data.trinaryReadsDF(
        df_sample=df)
    # See if number of "-1" is excessive
    dff = df_result + df_result.applymap(lambda v: -np.abs(v))
    frac_minus1 = -dff.sum().sum()  \
        /(2*len(df_result)*len(df_result.columns))
    self.assertLess(frac_minus1, 0.25)
    # Smoke tests for csv
    df_result = transform_data.trinaryReadsDF(
        csv_file="AM_MDM_Mtb_transcripts_DEseq.csv",
        is_display_errors=False,
        is_time_columns=False)

  # TODO: Fix so working with the same transformation of features,
  #       either all genes features or all gene-groups.
  def testTrinaryReadsDF2(self):
    return
    # Checks that trinary values computed directly from reads
    # are the same as those of normalized samples.
    # Get raw value of read counts
    provider = DataProvider()
    provider.do()
    #
    def calcTrinaryTimeSample(time_index):
        """
        Calculates the trinary value of a time sample
        :param str time_index: name of time value
        """
        int_index = int(time_index[1:])
        df0 = provider.dfs_read_count[0]
        num = len(provider.dfs_read_count)
        ser = pd.Series(np.repeat(0, len(df0.index)), index=df0.index)
        for idx in range(num):
            ser += provider.dfs_read_count[idx][int_index]
        df = pd.DataFrame(ser/num)
        df_result = transform_data.trinaryReadsDF(df_sample=df)
        return df_result.T
    #
    data = TrinaryData()
    data.df_X.columns = data.features
    for time_index in data.df_X.index:
      df_result = calcTrinaryTimeSample(time_index)
      import pdb; pdb.set_trace()
        
  def testCalcTrinaryComparison(self):
    if IGNORE_TEST:
      return
    df_in = pd.DataFrame({'a': [4, 0.20, 1]})
    df_expected = pd.DataFrame({'a': [1, -1, 0]})
    df_out = transform_data.calcTrinaryComparison(df_in)
    self.assertTrue(df_out.equals(df_expected))
    #
    df_out = transform_data.calcTrinaryComparison(df_in,
        ser_ref=df_in['a'])
    trues = [v == 0 for v in df_out['a']]
    self.assertTrue(all(trues))

  def testStripReplicaString(self):
    if IGNORE_TEST:
      return
    TIME = "TO"
    SIZE = 3
    names = ["%s.%d" % (TIME, n) for n in range(SIZE)]
    result = transform_data.stripReplicaString(names)
    self.assertEqual(result[0], TIME)
    self.assertEqual(len(result), SIZE)

  def testRemoveGenesWithExcessiveReplicationVariance(self):
    if IGNORE_TEST:
      return
    trinary = TrinaryData(is_averaged=False, is_dropT1=False,
        is_regulator=False)
    df_base = transform_data.removeGenesWithExcessiveReplicationVariance(
        trinary.df_X)
    for max_var in [1, 2, 3]:
      df = transform_data.removeGenesWithExcessiveReplicationVariance(
          trinary.df_X, max_var=max_var)
      self.assertGreaterEqual(len(df_base.columns), len(df.columns))

  def testConvertUnconvertToFromLog2(self):
    if IGNORE_TEST:
      return
    def test(pd_obj):
      if isinstance(pd_obj, pd.DataFrame):
        base_obj = DF
      else:
        base_obj = SER
      obj1 = transform_data.convertToLog2(base_obj)
      obj2 = transform_data.unconvertFromLog2(obj1)
      if isinstance(pd_obj, pd.DataFrame):
        ser2 = obj2["a"]
      else:
        ser2 = obj2
      ser2.loc[0] = 0
      trues = [np.isclose(v1, v2) for v1, v2 in zip(ser2, SER)]
      self.assertTrue(all(trues))
    #
    test(SER)
    test(DF)

    ser = transform_data.convertToLog2(SER)
    ser1 = transform_data.unconvertFromLog2(ser)
    ser1.loc[0] = 0
    trues = [np.isclose(v1, v2) for v1, v2 in zip(ser1, SER)]
    self.assertTrue(all(trues))
    
    

if __name__ == '__main__':
  unittest.main()
