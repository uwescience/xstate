"""Creates Data Used in Figures for Paper Submission"""

"""
TODO
1. makeCrossValidationData should be a tool that produces a CSV file.
2. Have a companion shell script so can run in parallel
3. Create a class so merge the result.
"""

from common import constants as cn
from common_python.util import dataframe
from tools.make_classification_data import ClassificationData

import argparse
import copy
import numpy as np
import os
import pandas as pd

# Constants
NUM_GENE = 11  # Number of genes considered
CV_CALCULATION_FILENAME = "cv_calculation"
DIR_PATH =  os.path.dirname(os.path.abspath(__file__))


class CrossValidationData(object):

  def __init__(self, num_clf=100, max_gene=NUM_GENE, dir_path=DIR_PATH):
    """
    Parameters
    ----------
    num_clf: int
        number of classifiers in ensemble
    max_gene: int
        Maximum number of genes considerd
    dir_path: str
        Directory where files are saved and read
    """
    self.num_clf = num_clf
    self.max_gene = max_gene
    self.dir_path = dir_path
    data = ClassificationData()
    self.namespace_dct = data.getDct()  # Serialized data from classifier constructions
    # Hypoxia data sources keyed by name of data source
    self.data_dct = self.namespace_dct["DATA_DCT"]
    # Collections of genes keyed by name of gene collection
    self.gene_dct = self.namespace_dct["GENE_DCT"]
    # Classifiers keyed by training data and genes
    self.classifier_dct = self.namespace_dct["CLASSIFIER_DCT"]

  def make(self, indices=None, num_iter=10):
    """
    Creates the data needed for accuracy plots based on cross validations.
  
    Parameters
    ----------
    indices: list-int
        Indicies of keys to process
    num_iter: int
        Number of iterations of cross validation
  
    Returns
    -------
    pd.DataFrame
        index: int
            maximum importance rank of the gene used to construct the classifier
        column: str (classifier name)
    """
    ranks = list(range(1, self.max_gene+1))
    columns = []
    result_dct = {}
    if indices is None:
      classifier_dct = CLASSIFIER_DCT
    else:
      keys = list(self.classifier_dct.keys())
      classifier_dct = {k: self.classifier_dct[k] for k in keys
          if keys.index(k) in indices}
    for key, clf in classifier_dct.items():
      classifier_name = "--".join(key)
      result_dct[classifier_name] = []
      columns.append(classifier_name)
      # Accuracy of classifier that includes this gene and all preceding genes
      accuracy_dct = {}
      trinary = copy.deepcopy(self.data_dct[key[0]])
      trinary.df_X = dataframe.subset(trinary.df_X, self.gene_dct[key[1]])
      for rank in ranks:
        accuracy = clf.crossValidate(
            trinary, num_iter=num_iter, num_holdout=1, filter_high_rank=rank,
                  size=self.num_clf)
        result_dct[classifier_name].append(accuracy)
    #
    df = pd.DataFrame(result_dct, index=ranks)
    df.index = [v + 1 for v in df.index]
    if self.dir_path is not None:
      path self._makeFilePath(indices)
      df.to_csv(path)
    return df

  def _makeFilePath(self, indices):
    """
    Constructs the path for the indices.

    Parameters
    ----------
    indices: list-int
    
    Returns
    -------
    str
    """
    sfx = "_".join([str(v) for v in indices])
    filename = "%s_%s.csv" % (CV_CALCULATION_FILENAME, sfx)
    return os.path.join(self.dir_path, filename)

  def _getFilePaths(self):
    """
    Gets the cross validation files in the directory.
 
    Returns
    -------
    list-str
    """
    files = os.listdir(self.dir_path)
    files = [CV_CALCULATION_FILENAME in f for f in files]
    return files

  def clean(self):
    """
    Removes all existing cross validation files.
    """
    ffiles = self._get_files()
    for ffile in ffiles:
      os.remove(ffile)

  def report(base_name):
    """
    Merges files with the same base name as columns in a dataframe.

    Parameters
    ----------
    base_name: str
    
    Returns
    -------
    pd.DataFrame
    """
    files = self._get_files()
    dfs = []
    for ffile in files:
      df = pd.csv_read(ffile)
      dfs.append(df)
    return pd.concat(dfs, axis=1)

if __name__ == '__main__':
  desc = 'Create cross validation data'
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('--start', '-s', help='starting index',
      default=-1, type=int)
  parser.add_argument('--end', '-e',  help='ending index',
      type=int, default=-1)
  parser.add_argument('--clean', '-c',  help='Clean the directory',
      type=int, default=7)
  args = parser.parse_args()
  #
  data = CrossValidationData()
  if args.clean:
    data.clean()
  if (start >= 0) and (end >= 0):    
    indices = list(range(args.start, args.end))
    _ = data.makeCrossValidationData(indices=indices)
