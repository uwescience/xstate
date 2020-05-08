import common.constants as cn
from common_python.testing import helpers
from common_python.util.persister import Persister
from classifier import  \
    main_multi_classifier_feature_optimizer as main

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_REPORT = False
FILENAME =  "test_main_multi_classifier.pcl"
DIR = os.path.dirname(os.path.abspath("__file__"))
PERSISTER_PATH = os.path.join(DIR, FILENAME)
MCFO_KWARGS = {
    "num_exclude_iter": 2,
    }


class TestFunctions(unittest.TestCase):

  def _remove(self):
    if os.path.exists(PERSISTER_PATH):
      os.remove(PERSISTER_PATH)

  def setUp(self):
    self._remove()

  def tearDown(self):
    self._remove()

  def testGetData(self):
    if IGNORE_TEST:
      return
    df_X, ser_y = main._getData()
    self.assertTrue(isinstance(df_X, pd.DataFrame))
    self.assertTrue(isinstance(ser_y, pd.Series))

  def testMakePath(self):
    if IGNORE_TEST:
      return
    path = main._makePath()
    try:
      fd = open(PERSISTER_PATH, "w")
      fd.close()
    except:
      self.assertTrue(False)

  def testRun(self):
    if IGNORE_TEST:
      return
    main.run(PERSISTER_PATH, True, max_iter=1,
        is_report=False, mcfo_kwargs=MCFO_KWARGS)
    persister = Persister(PERSISTER_PATH)
    self.assertTrue(persister.isExist())
    optimizer = persister.get()
    self.assertTrue(isinstance(
        optimizer.fit_result_dct, dict))
    #
    main.run(PERSISTER_PATH, False, max_iter=1,
        is_report=False, mcfo_kwargs=MCFO_KWARGS)
    optimizer2 = persister.get()
    for cls in optimizer.fit_result_dct.keys():
      self.assertTrue(
          len(optimizer.fit_result_dct[cls]) ==
          len(optimizer2.fit_result_dct[cls]))

  def testReport(self):
    if IGNORE_TEST:
      return
    main.run(PERSISTER_PATH, True, max_iter=1,
        is_report=False, mcfo_kwargs=MCFO_KWARGS)
    main.report(PERSISTER_PATH)

  def testWriteFitResultCSV(self):
    if IGNORE_TEST:
      return
    main.run(PERSISTER_PATH, True, max_iter=1,
        is_report=False, mcfo_kwargs=MCFO_KWARGS)
    df = main.makeFitResultCSV(path=PERSISTER_PATH,
        csv_path=None)
    self.assertTrue(helpers.isValidDataFrame(df,
        cn.FIT_RESULT_COLUMNS))


if __name__ == '__main__':
  unittest.main()
