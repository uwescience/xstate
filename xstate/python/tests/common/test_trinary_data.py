import common.constants as cn
from common.trinary_data import TrinaryData, NormalizedData
from common_python.testing import helpers

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
NUM_REPL = 3


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
      self.assertEqual(len(df_X), len(ser_y))
    #
    for cls in [TrinaryData, NormalizedData]:
      data1 = cls(is_averaged=False)
      data2 = cls(is_averaged=True)
      for data in [data1, data2]:
        test(data.df_X, data.ser_y)
      self.assertGreater(len(data1.df_X), len(data2.df_X))
      # Replicated data should have 3 Normoxia states
      self.assertEqual(data1.ser_y[
          data1.ser_y==0].count(),
          NUM_REPL)
    

if __name__ == '__main__':
  unittest.main()
