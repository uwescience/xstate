import common.constants as cn
from paper import cross_validation_data as cvd

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
DIR_PATH =  os.path.dirname(os.path.abspath(__file__))


class TestFunctions(unittest.TestCase):

  def setUp(self):
    self.data = cvd.CrossValidationData(num_clf=10, dir_path=DIR_PATH)

  def remove(self):
    paths = [f for f in os.listdir(DIR_PATH) if ("csv" in f) and ("cv_" in f)]
    for path in paths:
      os.remove(path)

  def testMakeData(self):
    if IGNORE_TEST:
      return
    indices = [0, 1]
    df = self.data.make(indices=indices, num_iter=2)
    self.assertTrue(isinstance(df, pd.DataFrame))
    self.assertGreater(len(df), 0)
    self.assertEqual(len(df.columns), len(indices))


if __name__ == '__main__':
  unittest.main()
