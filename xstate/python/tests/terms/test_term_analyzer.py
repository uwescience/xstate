from terms import term_analyzer
from common.data_provider import DataProvider
from common_python.testing import helpers
import common.constants as cn

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False


class TestTermAnalyzer(unittest.TestCase):

  def setUp(self):
    self.provider = DataProvider()
    self.provider.do()
    self.analyzer = term_analyzer.TermAnalyzer(
        self.provider.df_ec_terms, is_plot=IS_PLOT)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(helpers.isValidDataFrame(self.analyzer.df_term,
        self.analyzer.df_term.columns))

  def testPlotTermHeatmap(self):
    if IGNORE_TEST:
      return
    self.analyzer.plotTermHeatmap(is_plot=IS_PLOT)


if __name__ == '__main__':
  unittest.main()
