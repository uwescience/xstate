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


class TestDataTransformer(unittest.TestCase):

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
    df_result = transform_data.trinaryReadsDF(df_sample=df)
    # See if number of "-1" is excessive
    dff = df_result + df_result.applymap(lambda v: -np.abs(v))
    frac_minus1 = -dff.sum().sum()  \
        /(2*len(df_result)*len(df_result.columns))
    self.assertLess(frac_minus1, 0.25)
    # Smoke tests for csv
    df_result = transform_data.trinaryReadsDF(
        csv_file="AM_MDM_Mtb_transcripts_DEseq.csv",
        is_display_errors=False)

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


if __name__ == '__main__':
  unittest.main()
