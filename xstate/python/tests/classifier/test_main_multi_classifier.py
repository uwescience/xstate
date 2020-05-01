import common.constants as cn
from common_python.util.persister import Persister
from classifier import main_multi_classifier as main

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_REPORT = False
FILENAME =  "test_main_multi_classifier.pcl"
DIR = os.path.dirname(os.path.abspath("__file__"))
FILEPATH = os.path.join(DIR, FILENAME)


class TestFunctions(unittest.TestCase):

  def _remove(self):
    if os.path.exists(FILEPATH):
      os.remove(FILEPATH)

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
      fd = open(FILEPATH, "w")
      fd.close()
    except:
      self.assertTrue(False)

  def testRun(self):
    if IGNORE_TEST:
      return
    main.run(FILEPATH, True, max_iter=1,
        is_report=False)
    persister = Persister(FILEPATH)
    self.assertTrue(persister.isExist())
    optimizer = persister.get()
    self.assertTrue(isinstance(
        optimizer.feature_dct, dict))
    #
    main.run(FILEPATH, False, max_iter=1,
        is_report=False)
    optimizer2 = persister.get()
    for cls in optimizer.feature_dct.keys():
      self.assertTrue(
          optimizer.feature_dct[cls] ==
          optimizer2.feature_dct[cls])

  def testReport(self):
    if IGNORE_TEST:
      return
    main.run(FILEPATH, True, max_iter=1,
        is_report=False)
    main.report(FILEPATH)


if __name__ == '__main__':
  unittest.main()
