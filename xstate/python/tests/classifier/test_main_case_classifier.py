# TODO: Create CSV matrix files for existing data
from common_python.testing import helpers
from common_python.classifier import feature_analyzer
import classifier.main_case_classifier as main

import numpy as np
import os
import pandas as pd
import shutil
import unittest


IGNORE_TEST = False
IS_PLOT = False
DIR = os.path.dirname(os.path.abspath("__file__"))
TEST_OUT_PATH = os.path.join(DIR,
     "test_main_case_clasifier.csv")
SAMPLE_PATH = os.path.join(cn.DATA_DIR, "samples")
SAMPLE_PATH = os.path.join(SAMPLE_PATH,
    "galagan_raw_hypoxia_ts.csv")
TEST_IN_PATH = os.path.join(DIR,
     "test_main_case_clasifier.csv")


class TestFunctions(unittest.TestCase):

  def _remove(self):
    for path in [TEST_OUT_PATH]:
      if os.path.isfile(path):
        os.remove(path)

  def setUp(self):
    self._remove()

  def tearDown(self):
    self._remove()

  def testRun(self):
    if IGNORE_TEST:
      return
    return
    #
    with open(SAMPLE_PATH, "r") as fd:
      main.run(fd, TEST_OUT_PATH)
    self.assertTrue(os.path.isfile(TEST_OUT_PATH))


if __name__ == '__main__':
  unittest.main()
