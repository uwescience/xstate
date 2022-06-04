"""Creates Data Used in Figures for Paper Submission"""

"""
TODO
1. makeAccuracyData should be a tool that produces a CSV file.
2. Have a companion shell script so can run in parallel
3. Create a class so merge the result.
"""

from common import constants as cn
from common_python.util import dataframe
from tools.make_classification_data import ClassificationData

import copy
import numpy as np
import pandas as pd

# Constants
data = ClassificationData()
namespace_dct = data.getDct()  # Serialized data from classifier constructions
DATA_DCT = namespace_dct["DATA_DCT"]
GENE_DCT = namespace_dct["GENE_DCT"]
CLASSIFIER_DCT = namespace_dct["CLASSIFIER_DCT"]
NUM_GENE = 11  # Number of genes considered


def makeAccuracyData(num_clf=100, max_gene=NUM_GENE, path=None, indices=None):
  """
  Creates the data needed for accuracy plots.

  Parameters
  ----------
  num_clf: int
      number of classifiers in ensemble
  max_gene: int
      Maximum number of genes considerd
  path: str
      Save result to path
  indices: list-int
      Indicies of keys to process

  Returns
  -------
  pd.DataFrame
      index: int
          maximum importance rank of the gene used to construct the classifier
      column: str (classifier name)
  """
  ranks = list(range(1, max_gene+1))
  columns = []
  result_dct = {}
  if indices is None:
    classifier_dct = CLASSIFIER_DCT
  else:
    keys = list(CLASSIFIER_DCT.keys())
    classifier_dct = {k: CLASSIFIER_DCT[k] for k in keys
        if keys.index(k) in indices}
  for key, clf in classifier_dct.items():
    classifier_name = "--".join(key)
    result_dct[classifier_name] = []
    columns.append(classifier_name)
    # Accuracy of classifier that includes this gene and all preceding genes
    accuracy_dct = {}
    trinary = copy.deepcopy(DATA_DCT[key[0]])
    trinary.df_X = dataframe.subset(trinary.df_X, GENE_DCT[key[1]])
    for rank in ranks:
      accuracy = clf.crossValidate(
          trinary, num_iter=10, num_holdout=1, filter_high_rank=rank, size=num_clf)
      result_dct[classifier_name].append(accuracy)
  #
  df = pd.DataFrame(result_dct, index=ranks)
  df.index = [v + 1 for v in df.index]
  if path is not None:
    df.to_csv(path)
  return df
