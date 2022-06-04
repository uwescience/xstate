import common.constants as cn
from paper import make_data as md

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False


class TestFunctions(unittest.TestCase):

  def testMakeData(self):
    df = md.makeAccuracyData(indices=[0, 1])
    import pdb; pdb.set_trace()


if __name__ == '__main__':
  unittest.main()
