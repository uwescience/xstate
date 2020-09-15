import common.constants as cn
from common.data_provider import DataProvider
from common_python.util.persister import Persister
from tools import shared_data

import os
import numpy as np
import pandas as pd
import unittest

IGNORE_TEST = False
IS_PLOT = False
DIR = os.path.dirname(os.path.abspath(__file__))
PERSISTER_PATH = os.path.join(DIR,
    "persister_shared_data.pcl")


class TestSharedData(unittest.TestCase):

  def deleteFiles(self):
    if os.path.isfile(PERSISTER_PATH):
      os.remove(PERSISTER_PATH)

  def setUp(self):
    self.deleteFiles()
    self.persister = Persister(PERSISTER_PATH)

  def tearDown(self):
    self.deleteFiles()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    def test():
      data = shared_data.SharedData(persister=self.persister)
      self.assertTrue(isinstance(data.provider,
          DataProvider))
      self.assertTrue(isinstance(data.df_X, pd.DataFrame))
      self.assertTrue(isinstance(data.ser_y, pd.Series))
      self.assertTrue(isinstance(data.states,
          np.ndarray))
      self.assertEqual(len(data.states),
          len(data.collection_dct.keys()))
    # Test without persister
    test()
    self.assertTrue(self.persister.isExist())
    # Test with persister
    test()
      

if __name__ == '__main__':
  unittest.main()
