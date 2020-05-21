import common.constants as cn
from common_python.testing import helpers
from common_python.util.persister import Persister
from classifier import  \
    main_feature_equivalence_calculator as main
from common_python.classifier import  \
   multi_classifier_feature_optimizer as mcfo 

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_REPORT = False
PERSISTER_FILE =  \
     "test_main_feature_equivalence_calculator.pcl"
DIR = os.path.dirname(os.path.abspath("__file__"))
PERSISTER_PATH = os.path.join(DIR, PERSISTER_FILE)
OUT_FILE =  \
     "test_main_feature_equivalence_calculator.csv"
OUT_PATH = os.path.join(DIR, OUT_FILE)
FIT_RESULT_PATH = os.path.join(cn.DATA_DIR,
    "fit_result_tf.xlsx")
TEST_FIT_RESULT_PATH = os.path.join(cn.DATA_DIR,
    "fit_result_tf.xlsx")
STATE = 1


class TestFunctions(unittest.TestCase):

  def _remove(self):
    for path in [PERSISTER_PATH, OUT_PATH]:
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
    def test():
      persister = Persister(PERSISTER_PATH)
      self.assertTrue(persister.isExist())
      calculator = persister.get()
      self.assertTrue(isinstance(
          calculator.df_ria, pd.DataFrame))
      self.assertTrue(os.path.isfile(OUT_PATH))
    #
    main.run(STATE, 
        persister_path=PERSISTER_PATH, 
        out_path=OUT_FILE,
        is_restart=True,
        is_report=IS_REPORT,
        num_cross_iter=10,
        min_score=0.94)
    test()
    #
    main.run(STATE, 
        persister_path=PERSISTER_PATH, 
        out_path=OUT_FILE,
        is_restart=False,
        is_report=IS_REPORT,
        num_cross_iter=10,
        min_score=0.94)
    test()

  def testMakeFitResults(self):
    if IGNORE_TEST:
      return
    fit_results = \
       main.makeFitResult(1,
       TEST_FIT_RESULT_PATH, min_score=0.9)
    self.assertGreater(len(fit_results), 0)
    self.assertTrue(isinstance(fit_results[0],
        mcfo.FitResult))
    #
    fit_results = \
       main.makeFitResult(STATE,
       TEST_FIT_RESULT_PATH, min_score=1.0)
    self.assertEqual(len(fit_results), 0)


if __name__ == '__main__':
  unittest.main()
