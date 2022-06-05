import common.constants as cn
from paper.cross_validation_data import CrossValidationData

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False


class TestFunctions(unittest.TestCase):

  def setUp(self):
    self.data = CrossValidationData(num_clf=10)

  def testMakeData(self):
    df = self.data.make(indices=[0, 1], num_iter=2)
    import pdb; pdb.set_trace()


if __name__ == '__main__':
  unittest.main()
