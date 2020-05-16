import common.constants as cn
from common_python.testing import helpers
from common_python.util.persister import Persister
from classifier import  \
    main_feature_equivalence_calculator as main
from common_python.stateclassifier import  \
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
OUT_PATH = os.path.join(DIR, PERSISTER_FILE)
FIT_RESULT_PATH = os.path.join(cn.DATA_DIR,
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
    self.assertEqual(len(ser_y.values.unique()), 2)

  def testRun(self):
    if IGNORE_TEST:
      return
    fit_results =  \
        mcfo.MultiClassifierFeatureOptimizer.makeFitResult(
        FIT_RESULT_PATH,
        lambda r: r["idx"] == STATE)
    import pdb; pdb.set_trace()
    fit_results = fit_result[3]
    main.run(STATE, 
        persister_path=PERSISTER_PATH, 
        is_restart=True,
        is_report=False,
        fit_results=fit_results,
        num_cross_iter=2)
    persister = Persister(PERSISTER_PATH)
    self.assertTrue(persister.isExist())
    calculator = persister.get()
    self.assertTrue(isinstance(
        optimizer.ria_dct, dict))
    self.assertTrue(os.path.isfile(OUT_PATH))
    import pdb; pdb.set_trace()


if __name__ == '__main__':
  unittest.main()
