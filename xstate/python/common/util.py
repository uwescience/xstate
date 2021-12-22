""" Base utilities """

from common import constants as cn

import pandas as pd
import numpy as np


def convertToLog2(obj):
  """
  Converts to log2 units.

  Parameters
  ----------
  obj: DataFrame / Series
  
  Returns
  -------
  DataFrame
  """
  if isinstance(obj, pd.DataFrame):
    return obj.applymap(lambda v: np.log2(v)
        if v > cn.MIN_VALUE else np.log2(cn.MIN_VALUE))
  else:  # Series
    return obj.apply(lambda v: np.log2(v)
        if v > cn.MIN_VALUE else np.log2(cn.MIN_VALUE))

def unconvertFromLog2(obj):
  """
  Converts from log2 units.

  Parameters
  ----------
  obj: DataFrame / Series
  
  Returns
  -------
  DataFrame
  """
  if isinstance(obj, pd.DataFrame):
    return obj.applymap(lambda v: 2**v)
  else:  # Series
    return obj.apply(lambda v: 2**v)
