#!/usr/bin/env python
# coding: utf-8
"""
Constructs data used in the analysis of lab samples.
The data are saved in the Persister 
file $PROJECT_DIR/data/make_classification_data.pcl

Usage 1: Create the persister file
  data = ClassificationData()
  data.initialize()
  data.serialize()

Usage 2: Use an existing persister file
  data = ClassificationData()
  data.get(globals())
"""


from common import constants as cn
from common.trinary_data import TrinaryData, REF_TYPE_BIOREACTOR, \
    REF_TYPE_SELF, REF_TYPE_POOLED
from common.data_provider import DataProvider
from common import trinary_data
from common_python.plots import util_plots
from common_python.classifier import classifier_ensemble
from common_python.classifier import classifier_collection
from common_python.util.persister import Persister
from common_python.util import dataframe
from common import transform_data

import collections
import copy
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.decomposition import PCA
from sklearn import svm
from sklearn.model_selection import cross_val_score
import seaborn as sns

# Constants
DATA_FILE = "make_classification_data.pcl"
NUM_CLASSIFIER_IN_ENSEMBLE = 100



# Internal functions
def _updateSampleDct(sample_dct):
  """
  Transforms the elements of the sample dictionary
  into dictionaries.

  Parameters
  ----------
  sample_dct: dict
     key: ref_type (str) [.e.g, T0, POOLED]
     value: SampleData
  
  Returns
  -------
  dict:
     key: ref_type (str)
     value: dict
         key: sample name
         value: DataFrame
  """
  for name, sample_data in sample_dct.items():
    data_dct = {n: sample_data.__getattribute__(n)
        for n in trinary_data.SAMPLES}
    sample_dct[name] = data_dct


class ClassificationData():
  # Data preparation constants

  def __init__(self, persister_path=None):
    """
    Parameters
    ----------
    persister_path: str
        path to persister file
    """
    if persister_path is None:
      persister_path = os.path.join(cn.DATA_DIR, DATA_FILE)
    self.persister = Persister(persister_path)
    self.namespace_dct = {}  # Items that go in the caller's namespace

  def initialize(self):
    """
    Initializes the data. Defines and initializes all names added to globals().
    """
    #
    T0 = "T0"
    POOLED = "pooled"
    self._addName("T0", "T0")
    self._addName("POOLED", "pooled")
    self._addName("REF_TYPE_POOLED", REF_TYPE_POOLED)
    self._addName("REF_TYPE_BIOREACTOR", REF_TYPE_BIOREACTOR)
    self._addName("REF_TYPE_SELF", REF_TYPE_SELF)
    # Provider
    PROVIDER = DataProvider()
    self._addName("PROVIDER", PROVIDER)
    PROVIDER.do()
    TRINARY = TrinaryData()
    self._addName("TRINARY", TRINARY)
    # Gene Classes
    ALL_GENES = list(TRINARY.df_X.columns)
    self._addName("ALL_GENES", ALL_GENES)
    # Gene groupings
    MYCOBACTIN_GENES = [
      "Rv2377c",
      "Rv2378c",
      "Rv2379c",
      "Rv2380c",
      "Rv2381c",
      "Rv2382c",
      "Rv2383c",
      "Rv2384",
      "Rv2385",
      "Rv2386c",
    ]
    self._addName("MYCOBACTIN_GENES", MYCOBACTIN_GENES)
    BACTERIOFERRITIN_GENES = [
      "Rv2341", "Rv3841", 
    ]
    self._addName("BACTERIOFERRITIN_GENES", BACTERIOFERRITIN_GENES)
    MYCOBACTIN_BACTERIOFERRIN_GENES = list(MYCOBACTIN_GENES)
    self._addName("MYCOBACTIN_BACTERIOFERRIN_GENES", MYCOBACTIN_BACTERIOFERRIN_GENES)
    MYCOBACTIN_BACTERIOFERRIN_GENES.extend(BACTERIOFERRITIN_GENES)
    MYCOBACTIN_BACTERIOFERRITIN = "mycobactin_bacterioferritin"
    BACTERIOFERRITIN = "bacterioferritin"
    MYCOBACTIN = "mycobactin"
    ALL = "all"
    GENE_DCT = {MYCOBACTIN: MYCOBACTIN_GENES,
                BACTERIOFERRITIN: BACTERIOFERRITIN_GENES,
                MYCOBACTIN_BACTERIOFERRITIN: MYCOBACTIN_BACTERIOFERRIN_GENES,
                ALL: ALL_GENES,
               }
    GENE_GROUPS = list(GENE_DCT.keys())
    self._addName("GENE_GROUPS", GENE_GROUPS)
    for name in GENE_GROUPS:
      self._addName(name.upper(), name)
    self._addName("GENE_DCT", GENE_DCT)
    # Define the stage names
    STAGE_NAMES = list(cn.STATE_NAMES)
    self._addName("STAGE_NAMES", STAGE_NAMES)
    STAGE_NAMES.remove("Normoxia")
    STAGE_NAMES = np.array(STAGE_NAMES)
    # Bioreactor data calculated with two different references
    DATA_DCT = {
        T0: TrinaryData(is_regulator=False, is_dropT1=True,
                        is_averaged=True),
        POOLED: TrinaryData(is_regulator=False, is_dropT1=True,
                            is_averaged=True, 
                            calcRef=PROVIDER.calcRefPooled)
    }
    self._addName("DATA_DCT", DATA_DCT)
    SER_Y_DCT = {k: t.ser_y for k,t in DATA_DCT.items()}
    self._addName("SER_Y_DCT", SER_Y_DCT)
    # Feature vectors are specific to the gene subsets
    DF_X_DCT = {k: t.df_X.copy() for k,t in DATA_DCT.items()}
    DF_X_DCT = {k: df[MYCOBACTIN_GENES] for k, df in DF_X_DCT.items()}
    self._addName("DF_X_DCT", DF_X_DCT)
    # Sample data
    SAMPLE_DCT = {r: trinary_data.getSampleData(ref_type=r, is_regulator=False)
        for r in [REF_TYPE_BIOREACTOR, REF_TYPE_SELF, REF_TYPE_POOLED]}
    self._addName("SAMPLE_DCT", SAMPLE_DCT)
    # Classifiers
    max_rank = len(MYCOBACTIN_BACTERIOFERRIN_GENES)
    CLASSIFIER_BASE = classifier_ensemble.ClassifierEnsemble(
          classifier_ensemble.ClassifierDescriptorSVM(),
          filter_high_rank=max_rank, size=NUM_CLASSIFIER_IN_ENSEMBLE)
    self._addName("CLASSIFIER_BASE", CLASSIFIER_BASE)
    CLASSIFIER_DCT = {}
    self._addName("CLASSIFIER_DCT", CLASSIFIER_DCT)
    for trinary_key, trinary in DATA_DCT.items():
      for gene_key, gene_list in GENE_DCT.items():
        classifier = copy.deepcopy(CLASSIFIER_BASE)
        # Not all genes may be present in TrinaryData since they may be correlated or unvarying.
        df_X = dataframe.subset(trinary.df_X, gene_list, axis=1)
        classifier.fit(df_X, trinary.ser_y, class_names=STAGE_NAMES)
        CLASSIFIER_DCT[(trinary_key, gene_key)] = classifier
    _updateSampleDct(SAMPLE_DCT)
    # Construct derivative structures    
    self._addName("DF_X", DF_X_DCT[T0])
    self._addName("SER_Y", SER_Y_DCT[T0])
    self._addName("SAMPLE_DATA_DCT", SAMPLE_DCT[REF_TYPE_BIOREACTOR])
    self._addName("CLASSIFIER", CLASSIFIER_DCT[('T0', 'mycobactin')])
    key = (T0, "mycobactin_bacterioferritin")
    self._addName("GENES", CLASSIFIER_DCT[key].features)

  def _addName(self, name, value):
    """
    Adds the name and value to the namespace.
 
    Parameters
    ----------
    name: str
    value: object
    -------
    """
    stmt = "self.namespace_dct['%s'] = value" % name
    exec(stmt)

  def serialize(self):
    """
    Writes the current contents of self.namespace_dct to the persister.
    """
    self.persister.set(self.namespace_dct)

  def deserialize(self):
    """
    Recovers previously serialized data, initializing self.namespace_dct.
    -------
    """
    if not self.persister.isExist():
      raise ValueError("Persister file %s does not exist. Use serialize first."
          % self.persister.path)
    self.namespace_dct = self.persister.get()
    return self.namespace_dct

  def setNamespace(self, globals_dct):
    """
    Sets the globals provided based on the initialized namespace.

    Parameters
    ----------
    globals_dct: dict
    """
    for name, value in self.namespace_dct.items():
      globals_dct[name] = value

  def get(self, globals_dct):
    """
    Deserializes an existing persister file and initializes the namespace.

    Parameters
    ----------
    globals_dct: dict
    """
    self.deserialize()
    self.setNamespace(globals_dct)

      
if __name__ == '__main__':
  data = ClassificationData()
  data.initialize()
  data.serialize()
  print("Updated classification data in %s" % data.persister.path)
