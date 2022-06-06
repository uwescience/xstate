import common.constants as cn
from tools import cross_validation_data as cvd

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
DIR_PATH =  os.path.dirname(os.path.abspath(__file__))


class TestFunctions(unittest.TestCase):

  def setUp(self):
    self.remove()
    self.data = cvd.CrossValidationData(num_clf=10, dir_path=DIR_PATH)

  def tearDown(self):
    self.remove()

  def remove(self):
    paths = [os.path.join(DIR_PATH, f)
        for f in os.listdir(DIR_PATH) if ("csv" in f) and ("cv_" in f)]
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

  def testMakePath(self):
    if IGNORE_TEST:
      return
    path = self.data._makePath([0, 1])
    self.assertTrue("0_1" in str(path))
    path = self.data._makePath([2, 4])
    self.assertTrue("2_4" in str(path))

  def testGetPaths(self):
    if IGNORE_TEST:
      return
    path = self.data._makePath([0, 1])
    path.touch()
    path = self.data._makePath([2, 1])
    path.touch()
    paths = self.data._getPaths()
    self.assertEqual(len(paths), 2)

  def testClean(self):
    if IGNORE_TEST:
      return
    path = self.data._makePath([0, 1])
    path.touch()
    path = self.data._makePath([2, 1])
    path.touch()
    #
    self.data.clean()
    paths = self.data._getPaths()
    self.assertEqual(len(paths), 0)

  def test_dataframe(self):
    if IGNORE_TEST:
      return
    df0 = self.data.make(indices=[0], num_iter=2)
    df1 = self.data.make(indices=[1], num_iter=2)
    df = self.data.dataframe
    self.assertEqual(len(df.columns), 2)
    self.assertEqual(len(df), cvd.NUM_GENE)

  def test_dataframeDataDirectory(self):
    if IGNORE_TEST:
      return
    self.data = cvd.CrossValidationData()
    df = self.data.dataframe
    self.assertEqual(len(df), cvd.NUM_GENE)


if __name__ == '__main__':
  unittest.main()
