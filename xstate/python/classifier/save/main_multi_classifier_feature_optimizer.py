'''Runs MultiClassifier and reports results.'''
"""
MultiClassifierFeatureOptimizer is run with the
following hyperparameters
  num_iteration: number of genes
  feature_selector: default
  desired_accuracy: 0.01 (difference from using all features)
"""


import common.constants as cn
from common_python.classifier import  \
    multi_classifier_feature_optimizer as mcfo
from common_python.types.extended_dict  \
    import ExtendedDict
from common_python.classifier import feature_collection
from common_python.util.persister import Persister
from common.trinary_data import TrinaryData
from common.data_provider import DataProvider

import argparse
import os
import pandas as pd


PERSISTER_FILE =  \
    "main_multi_classifier_feature_optimizer.pcl"
FIT_RESULT_PATH = os.path.join(cn.DATA_DIR,
    "fit_result.csv")
MAX_ITER = 5000
MIN_INCR_SCORE = 0.01
DESIRED_ACCURACY = 1.0
NUM_HOLDOUTS = 1  # Holdouts in cross validation
NUM_CROSS_ITER = 150  # Cross validation iterations
BCFO_KWARGS = {
    "desired_accuracy": DESIRED_ACCURACY,
    "max_iter": MAX_ITER,
    "num_holdouts": NUM_HOLDOUTS,
    "num_cross_iter": NUM_CROSS_ITER,
    }
NUM_CLASSIFICATION_GROUPS = 30
IS_ONLY_TFS = True  # Only process transciption factors


def getFitResultFromPersister(path=None):
  """
  Returns the dictionary of FitResult by class.
  """
  persister = getPersister(path=path)
  optimizer = persister.get()
  return optimizer.fit_result_dct

def getPersister(path=None):
  if path is None:
    path = _makePath(filename=PERSISTER_FILE)
  return Persister(path)

def _makePath(filename=PERSISTER_FILE):
  this_dir = os.path.dirname(os.path.abspath(__file__))
  return os.path.join(this_dir, filename)

def _getData():
  provider = DataProvider()
  provider.do()
  trinary = TrinaryData(is_averaged=False,
      is_dropT1=False)
  if IS_ONLY_TFS:
    columns = set(trinary.df_X.columns).intersection(
        provider.tfs)
  else:
    columns = trinary.df_X.columns
  columns = list(columns)
  return trinary.df_X[columns], trinary.ser_y

def run(path, is_restart, max_iter=None,
    mcfo_kwargs={}, is_report=True):
  """
  Runs feature selection.
  :param str path: path to pcl file
  :param bool is_restart: Start a new analysis
  """ 
  df_X, ser_y = _getData()
  if max_iter is None:
    max_iter = len(df_X.columns) 
  persister = getPersister(path)
  if is_restart:
    if os.path.isfile(path):
      os.remove(path)
  if persister.isExist():
    if is_report:
      print ("\n***PCL found: %s\n" % path)
    optimizer = persister.get()
  else:
    if is_report:
      print ("\n***Creating new PCL file: %s\n"
          % path)
    bcfo_kwargs = dict(BCFO_KWARGS)
    bcfo_kwargs["max_iter"]  = max_iter
    optimizer = mcfo.MultiClassifierFeatureOptimizer(
        feature_collection_cl=  \
        feature_collection.FeatureCollection,
        persister=persister,
        **mcfo_kwargs, bcfo_kwargs=bcfo_kwargs)
  #
  optimizer.fit(df_X, ser_y)

def makeFitResultCSV(path=None,
    csv_path=FIT_RESULT_PATH):
  """
  Writes a CSV file of the fitness results to
  the data directory.
  :param str path: persister path
  :param str csv_path: path for writing CSV
      default: written to data directory
      None: no file written
  :return pd.DataFrame: Columns
      STATE: expression state
      GROUP: index for the sufficient feature set
      GENE_ID: Gene
      SCORE: cross validation score
      COUNT: Number of gene evaluations done
  """
  fit_result_dct = getFitResultFromPersister(
      path=path)
  result_dct = {k: [] for k in cn.FIT_RESULT_COLUMNS}
  for state in fit_result_dct.keys():
    for fit_result in fit_result_dct[state]:
      for gene in fit_result.sels:
        result_dct[cn.STATE].append(state)
        result_dct[cn.GROUP].append(
            fit_result.idx)
        result_dct[cn.GENE_ID].append(gene)
        result_dct[cn.SCORE].append(
            fit_result.sels_score)
        result_dct[cn.COUNT].append(
            fit_result.n_eval)
  df = pd.DataFrame(result_dct)
  if csv_path is not None:
    df.to_csv(csv_path, index=False)
  return df

def report_progress(path=None):
  """
  Reports progress on constructing sufficient feature sets.
  :param bool is_restart: Start a new analysis
  """
  def prt(header, dct):
    edct = ExtendedDict(dct)
    print("\n%s\n" % header)
    print(str(edct))
  #
  persister = getPersister(path)
  if persister.isExist():
    optimizer = persister.get()
    sels_len = [len(r) for 
        r in optimizer.fit_result_dct.values()]
    prt("\n**Sufficient feature sets constructed:\n",
        {s: len(optimizer.fit_result_dct[s])
        for s in optimizer.fit_result_dct.keys()})
    if False:
      if all([l > 0 for l in sels_len]):
        sels_dct = {c: r[0].sels
            for c, r in optimizer.fit_result_dct.items()}
        prt("\n**Selected features by state:\n", sels_dct)
        sels_score_dct = {c: r[0].sels_score
            for c, r in optimizer.fit_result_dct.items()}
        prt("\n**Scores for selected features by state:\n",
            sels_score_dct)
        all_score_dct = {c: r[0].all_score
            for c, r in optimizer.fit_result_dct.items()}
        prt("\n**Scores for all-features by state:\n",
            all_score_dct)
      else:
        print("***No data has been produced.")
  else:
    print("***Persister file not found: %s" % path)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      description="Run Feature Optimizer")
  parser.add_argument("-f", type=str,
      help="Pickle file with state",
      default=PERSISTER_FILE)
  parser.add_argument("--restart",
      action="store_true",
      help="Re-start the run from beginning",
      )
  parser.add_argument("--progress",
      action="store_true",
      help="Report run progress."
      )
  parser.add_argument("--write_csv",
      action="store_true",
      help="Write a CSV file of the fitness results."
      )
  args = parser.parse_args()
  path = _makePath(args.f)
  if args.progress:
    report_progress(path=path)
  elif args.write_csv:
    makeFitResultCSV()
  else:
    run(path=path,
        mcfo_kwargs=  \
        {"num_exclude_iter": 
        NUM_CLASSIFICATION_GROUPS},
        is_restart=args.restart)
    report_progress(path=path)
    makeFitResultCSV()
