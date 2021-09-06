import common.constants as cn
from common_python.testing import helpers
from common_python.classifier import feature_analyzer
import classifier.main_feature_analyzer as main

import numpy as np
import os
import pandas as pd
import shutil
import unittest


IGNORE_TEST = False
IS_PLOT = False
IS_REPORT = False
DIR = os.path.dirname(os.path.abspath("__file__"))
STATE = 1
# Create output paths with a variable for state
TEST_OUT_DIR_PAT = os.path.join(DIR,
    "test_main_feature_analyzer_%d")
FEATURE1 = "Rv1927"
FEATURE2 = "Rv1129c"
FEATURES = [FEATURE1, FEATURE2]


class TestFunctions(unittest.TestCase):

  def _remove(self):
    for metric in feature_analyzer.METRICS:
      path = TEST_OUT_DIR_PAT  % STATE
      if os.path.isdir(path):
        shutil.rmtree(path)

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
    main.run(STATE, out_dir_pat=TEST_OUT_DIR_PAT, is_restart=True,
        columns=FEATURES, num_cross_iter=2)
    for metric in feature_analyzer.METRICS:
      path = os.path.join(TEST_OUT_DIR_PAT % STATE,
          "%s.csv" % metric)
      self.assertTrue(os.path.isfile(path))

  def testRunStatus(self):
    if IGNORE_TEST:
      return
    main.run(STATE, out_dir_pat=TEST_OUT_DIR_PAT, is_restart=True,
        is_status=True, is_report=IS_PLOT,
        columns=FEATURES, num_cross_iter=2)


if __name__ == '__main__':
  unittest.main()
