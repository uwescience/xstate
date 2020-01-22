"""Uses a support vector machine (SVM) state classification."""

import common.constants as cn

from common.trinary_data import TrinaryData
from common_python.plots import util_plots

import numpy as np
import pandas as pd
from sklearn import svm
from sklearn.model_selection import cross_val_score


class SVMClassifier(object):

  def __init__(self, data=None):
    """
    :param NormalizedData data:
    """
    if data is None:
      self.data = TrinaryData()
    else:
      self.data = data
    self.scores = None

  def evaluate(self, iterations=5):
    """
    Does cross validations to assess accuracy of classifier.
    :param int iterations: Number of cross validations done
    Instance variables
      scores changed
    """
    indices = list(self.data.df_X.index)
    size = 5
    self.scores = []
    for _ in range(size):
      indices = np.random.permutation(indices)
      df_X = self.data.df_X.loc[indices, :]
      ser_y = self.data.ser_y.loc[indices]
      lin_clf.score(df_X, ser_y)
      self.scores.append(cross_val_score(lin_clf, df_X, ser_y, cv=3))
