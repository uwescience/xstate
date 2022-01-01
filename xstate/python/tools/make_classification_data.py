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
from common import sample_data
from common.sample_data import REF_TYPE_BIOREACTOR, \
    REF_TYPE_SELF, REF_TYPE_POOLED
from common.trinary_data import TrinaryData
from common.data_provider import DataProvider
from common.sample_data import SampleData
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
NUM_FEATURES = list(range(1, 13))  # Number of features in classifier plots
# Dataframe columns
COL_REF = "ref"
COL_GENE_GROUP = "Gene_GROUP"
COL_NUM_FEATURE = "num_feature"
COL_MEAN_ACCURACY = "mean_accuracy"
COL_STD_ACCURACY = "std_accuracy"
DF_ACCURACY_COLUMNS = [COL_REF, COL_GENE_GROUP, COL_NUM_FEATURE,
    COL_MEAN_ACCURACY, COL_STD_ACCURACY]


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
    SAMPLE_DCT = {r: sample_data.getSampleData(ref_type=r, is_regulator=False)
        for r in [REF_TYPE_BIOREACTOR, REF_TYPE_SELF, REF_TYPE_POOLED]}
    self._addName("SAMPLE_DCT", SAMPLE_DCT)
    SAMPLE_AVG_DCT = {r: sample_data.getSampleData(ref_type=r,
        is_regulator=False, is_average=True)
        for r in [REF_TYPE_BIOREACTOR, REF_TYPE_SELF, REF_TYPE_POOLED]}
    self._addName("SAMPLE_AVG_DCT", SAMPLE_AVG_DCT)
    # Classifiers
    num_feature = len(MYCOBACTIN_BACTERIOFERRIN_GENES)
    CLASSIFIER_BASE = classifier_ensemble.ClassifierEnsemble(
          classifier_ensemble.ClassifierDescriptorSVM(),
          filter_high_rank=num_feature, size=NUM_CLASSIFIER_IN_ENSEMBLE)
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
    # Construct derivative structures    
    self._addName("DF_X", DF_X_DCT[T0])
    self._addName("SER_Y", SER_Y_DCT[T0])
    self._addName("SAMPLE_DATA_DCT", SAMPLE_DCT[REF_TYPE_BIOREACTOR])
    self._addName("CLASSIFIER", CLASSIFIER_DCT[('T0', 'mycobactin')])
    key = (T0, "mycobactin_bacterioferritin")
    self._addName("GENES", CLASSIFIER_DCT[key].features)
    # Accuracy calculations for classifiers
    DF_ACCURACY = self.calcAccuracy()
    self._addName("DF_ACCURACY", DF_ACCURACY)

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

  def calcAccuracy(self, num_features=NUM_FEATURES, num_clf=100, is_debug=False):
    """
    Calculates the accuracy of classifiers using 10 iterations of 
    cross validation with one holdout per state (stage).

    Parameters
    ----------
    num_features: list-int
    num_clf: number of classifiers in the ensemble
    is_debug: bool
        Creates dummy data 
    
    Returns
    -------
    DataFrame:
        COL_REF: how reference is calculated for gene expressions
        COL_GENE_GROUP: grouping of genes used in classifier
        COL_NUM_FEATURE: number of features in classifiers
        COL_MEAN_ACCURACY: mean accuracy of the classifiers
        COL_STD_ACCURACY: standard deviation of accuracy
    """
    classifier_dct = self.namespace_dct["CLASSIFIER_DCT"]
    data_dct = self.namespace_dct["DATA_DCT"]
    gene_dct = self.namespace_dct["GENE_DCT"]
    line_dct = {r: l for r, l in zip(data_dct.keys(), ["-", "--"])}
    accuracy_dct = {c: [] for c in DF_ACCURACY_COLUMNS}
    for (ref, group), clf in classifier_dct.items():
      num_features = list(range(1, 13)) 
      num_features.insert(0, 1)
      trinary = copy.deepcopy(data_dct[ref])
      trinary.df_X = dataframe.subset(trinary.df_X, gene_dct[group])
      for num_feature in num_features:
        if is_debug:
          # Create a dummy value
          mean_accuracy = np.random.rand()
        else:
          mean_accuracy = clf.crossValidate(
              trinary, num_iter=10, num_holdout=1,
              filter_high_rank=num_feature, size=num_clf)
        accuracy_dct[COL_REF].append(ref)
        accuracy_dct[COL_GENE_GROUP].append(group)
        accuracy_dct[COL_NUM_FEATURE].append(num_feature)
        accuracy_dct[COL_MEAN_ACCURACY].append(mean_accuracy)
        std_accuracy = np.sqrt(mean_accuracy*(1 - mean_accuracy)/num_clf)
        accuracy_dct[COL_STD_ACCURACY].append(std_accuracy)
    df_accuracy = pd.DataFrame(accuracy_dct)
    return df_accuracy

      
if __name__ == '__main__':
  data = ClassificationData()
  data.initialize()
  data.serialize()
  print("Updated classification data in %s" % data.persister.path)
