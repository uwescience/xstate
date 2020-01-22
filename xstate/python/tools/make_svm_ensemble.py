"""Creates an SVM classifier from the time course expression data."""

import add_path
from common import constants as cn
from common.trinary_data import TrinaryData
from common_python.classifier import classifier_ensemble
from common_python.util import persister
from common import transform_data

import argparse
import os
import numpy as np
import pandas as pd

import numpy as np
import pandas as pd

NUM_FEATURES = 15
NUM_CLASSIFIERS = 50

def make(num_features=NUM_FEATURES,
    num_classifiers=NUM_CLASSIFIERS, 
    file_path=cn.ENSEMBLE_PATH, is_force=True):
  """
  Creates a classifier ensemble from the time course data,
  and fits the classifier.
  :param int num_features: number of features in the ensemble
  :param int num_classifiers: number of classifiers
  :param str file_path: path where classifier is exported
  :param bool is_force: Create new classifier even
      if one exists already.
  :return  ClassifierEnsemble:
  """
  # See if there's a classifier already
  if not is_force:
    try:
      return classifier_ensemble.ClassifierEnsemble.deserialize(
          file_path)
    except ValueError:
      pass
  # Construct a new classifier
  svm_ensemble = classifier_ensemble.ClassifierEnsemble(
          classifier_ensemble.ClassifierDescriptorSVM(), 
          filter_high_rank=num_features, size=num_classifiers)
  data = TrinaryData()
  data.df_X.columns = data.features
  svm_ensemble.fit(data.df_X, data.ser_y)
  svm_ensemble.serialize(file_path)
  return svm_ensemble

if __name__ == '__main__':
  # Do arg parse with errors
  desc = 'Create classifier ensemble in a destination file'
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('--dest', '-d', help='destination file',
      default=cn.ENSEMBLE_PATH)
  parser.add_argument('--numf', '-f', 
      help='number of features in classifiers',
      type=int, default=NUM_FEATURES)
  parser.add_argument('--numc', '-c', 
      help='number of classifiers in ensemble',
      type=int, default=NUM_CLASSIFIERS)
  args = parser.parse_args()
  make(num_features=args.numf, num_classifiers=args.numc,
      file_path=args.dest)
  print("Success!")
