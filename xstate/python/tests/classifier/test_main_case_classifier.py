# TODO: Create CSV matrix files for existing data
import common.constants as cn
import common_python.constants as ccn
from common_python.testing import helpers
from common_python.classifier import feature_analyzer
import classifier.main_case_classifier as main

import numpy as np
import os
import pandas as pd
import shutil
import unittest


IGNORE_TEST = True
IS_PLOT = True
DIR = os.path.dirname(os.path.abspath(__file__))
TEST_OUT_PATH = os.path.join(DIR,
     "test_main_case_clasifier.csv")
TEST_IN_PATH = os.path.join(cn.TRINARY_SAMPLES_DIR,
     "sherman.csv")
STATE = 1


class TestFunctions(unittest.TestCase):

  def _remove(self):
    for path in [TEST_OUT_PATH]:
      if os.path.isfile(path):
        os.remove(path)

  def setUp(self):
    self._remove()

  def tearDown(self):
    self._remove()

  def testRunState(self):
    if IGNORE_TEST:
      return
    df_instance = pd.read_csv(TEST_IN_PATH)
    arguments = main.Arguments(
        state=STATE, df=df_instance, num_fset=10)
    df = main._runState(arguments)
    columns = expected_columns=[ccn.FEATURE_VECTOR,
        ccn.SIGLVL, cn.STATE, main.INSTANCE]
    self.assertTrue(helpers.isValidDataFrame(df,
        expected_columns=columns,
        nan_columns=columns))
  
  def testRun(self):
    # TESTING
    #
    with open(TEST_IN_PATH, "r") as fd:
      main.run(fd, TEST_OUT_PATH, num_fset=2)
    self.assertTrue(os.path.isfile(TEST_OUT_PATH))
    self.assertTrue(os.path.isfile(TEST_OUT_PATH))


if __name__ == '__main__':
  unittest.main()
