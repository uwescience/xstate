from common import util
import common.constants as cn
from common_python.testing import helpers

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
SER = pd.Series([float(v) for v in range(10)])
DF = pd.DataFrame({"a": SER})
DF["b"] = DF["a"]*10


class TestFunctions(unittest.TestCase):

  def testConvertUnconvertToFromLog2(self):
    if IGNORE_TEST:
      return
    def test(pd_obj):
      if isinstance(pd_obj, pd.DataFrame):
        base_obj = DF
      else:
        base_obj = SER
      obj1 = util.convertToLog2(base_obj)
      obj2 = util.unconvertFromLog2(obj1)
      if isinstance(pd_obj, pd.DataFrame):
        ser2 = obj2["a"]
      else:
        ser2 = obj2
      ser2.loc[0] = 0
      trues = [np.isclose(v1, v2) for v1, v2 in zip(ser2, SER)]
      self.assertTrue(all(trues))
    #
    test(SER)
    test(DF)

    ser = util.convertToLog2(SER)
    ser1 = util.unconvertFromLog2(ser)
    ser1.loc[0] = 0
    trues = [np.isclose(v1, v2) for v1, v2 in zip(ser1, SER)]
    self.assertTrue(all(trues))
    
    

if __name__ == '__main__':
  unittest.main()
