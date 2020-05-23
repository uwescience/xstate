import common.constants as cn
from common_python.testing import helpers
import classifier.main_feature_analyzer as main

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_REPORT = False
DIR = os.path.dirname(os.path.abspath("__file__"))
STATE = 1
TEST_OUT_PATH_PAT = {}
# Create output paths with a variable for state
TEST_OUT_PATH_PAT = os.path.join(DIR,
    "test_main_feature_analyzer_%s_%d.csv")
FEATURE1 = "Rv0158"
FEATURE2 = "Rv1460"
FEATURES = [FEATURE1, FEATURE2]


class TestFunctions(unittest.TestCase):

  def _remove(self):
    for metric in main.METRICS:
      path = TEST_OUT_PATH_PAT  % (metric, STATE)
      if os.path.exists(path):
        os.remove(path)

  def setUp(self):
    self._remove()

  def tearDown(self):
    self._remove()

  def testGetData(self):
    if IGNORE_TEST:
      return
    df_X, ser_y = main._getData(STATE)
    self.assertTrue(isinstance(df_X, pd.DataFrame))
    self.assertTrue(isinstance(ser_y, pd.Series))
    self.assertEqual(len(ser_y.unique()), 2)

  def testRun(self):
    if IGNORE_TEST:
      return
    #
    main.run(STATE, out_path_pat=TEST_OUT_PATH_PAT,
        is_report=IS_REPORT,
        columns=FEATURES,
        num_cross_iter=2)
    for metric in main.METRICS:
      self.assertTrue(os.path.isfile(
          TEST_OUT_PATH_PAT % (metric, STATE)))


if __name__ == '__main__':
  unittest.main()
