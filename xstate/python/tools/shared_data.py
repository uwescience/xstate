"""Common Data Used in Tools."""

import common.constants as cn
from common.data_provider import DataProvider
import common_python.classifier.feature_analyzer as fa
from common_python.classifier  \
    import feature_set_collection as fsc
from common_python.classifier.feature_set import FeatureSet
from common import trinary_data
from common_python.util.persister import Persister

import argparse
import numpy as np
import os
import pandas as pd

DIR = os.path.dirname(os.path.abspath(__file__))
PERSISTER_PATH = os.path.join(DIR,
    "persister_shared_data.pcl")
PERSISTER = Persister(DIR)


class SharedData(object):

  def __init__(self, persister=PERSISTER):
    if persister.isExist():
      shared_data = persister.get()
      for key in shared_data.__dict__.keys():
        self.__setattr__(
            key, shared_data.__getattribute__(key))
    else:
      self.provider = DataProvider()
      self.trinary = trinary_data.TrinaryData(
          is_averaged=False, is_dropT1=False,
          is_regulator=True)
      self.df_X = self.trinary.df_X
      self.ser_y = self.trinary.ser_y
      self.states = self.ser_y.unique()
      data_dir = os.path.join(cn.DATA_DIR,
           "feature_analyzer")
      data_path_pat = os.path.join(data_dir, "%d") 
      data_dir = os.path.join(cn.DATA_DIR,
           "feature_analyzer")
      analyzer_dct = fa.deserialize({s: data_path_pat % s for s in self.states})
      analyzers = analyzer_dct.values()
      self.collection_dct = {s: 
          fsc.FeatureSetCollection.deserialize(
          data_path_pat % s) for s in self.states}
      _ = [c.df_fv for c in self.collection_dct.values()]
      persister.set(self)
