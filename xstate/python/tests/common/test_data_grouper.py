from common.data_grouper import DataGrouper
import common.constants as cn
from common.data_provider import DataProvider
from common_python.testing import helpers

import copy
import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False


class TestFunction(unittest.TestCase):

  def setUp(self):
    self.grouper = DataGrouper()

  def testGroupTrinaryValueRows(self):
    if IGNORE_TEST:
      return
    self.grouper._makeGroupDF()
    expecteds = set([cn.CNT_GROUP, cn.CNT_UP, cn.CNT_DOWN,
        cn.GENE_IDS])
    difference = expecteds.symmetric_difference(self.grouper.df_group)
    helpers.isValidDataFrame(self.grouper.df_group,
        [cn.CNT_GROUP, cn.CNT_UP, cn.CNT_DOWN, cn.GENE_IDS])

  def testPruneSize(self):
    if IGNORE_TEST:
      return
    self.grouper._makeGroupDF()
    length = len(self.grouper.df_group)
    self.grouper._pruneSize(min_size=0)
    self.assertEqual(length, len(self.grouper.df_group) + 1)
    self.grouper._pruneSize(min_size=1)
    self.assertLess(len(self.grouper.df_group), length)

  def testMakeGeneGroupDF(self):
    if IGNORE_TEST:
      return
    self.grouper._makeGeneGroupDF()
    helpers.isValidDataFrame(self.grouper.df_gene_group,
        [cn.GROUP])
    df = copy.deepcopy(self.grouper.df_gene_group)
    self.grouper._makeGeneGroupDF(is_include_0_group=True)
    helpers.isValidDataFrame(self.grouper.df_gene_group,
        [cn.GROUP])
    self.assertGreater(len(self.grouper.df_gene_group), len(df))
    
    

  def testDo(self):
    # Smoke test
    if IGNORE_TEST:
      return
    self.grouper.do()
    

if __name__ == '__main__':
  unittest.main()
