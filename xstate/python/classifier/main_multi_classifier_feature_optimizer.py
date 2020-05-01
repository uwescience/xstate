'''Runs MultiClassifier and reports results.'''
"""
MultiClassifier is run with the following hyperparameters
  num_iteration: number of genes
  feature_selector: default
  max_degrade: 0.01 (difference from using all features)
"""


from common_python.classifier import  \
    multi_classifier_feature_optimizer as mcfo
from common_python.types.extended_dict  \
    import ExtendedDict
from common_python.classifier import feature_collection
from common_python.util.persister import Persister
from common.trinary_data import TrinaryData

import argparse
import os


SERIALIZE_FILE =  \
    "main_multi_classifier_feature_optimizer.pcl"
MAX_ITER = 2500
MIN_INCR_SCORE = 0.01
MAX_DEGRADE = 0.01
NUM_HOLDOUTS = 1  # Holdouts in cross validation
NUM_CROSS_ITER = 20  # Cross validation iterations
BCFO_KWARGS = {
    "max_degrade": MAX_DEGRADE,
    "max_iter": MAX_ITER,
    "num_holdouts": NUM_HOLDOUTS,
    "num_cross_iter": NUM_CROSS_ITER,
    }


def _getPersister(path=None):
  if path is None:
    path = _makePath(filename=SERIALIZE_FILE)
  return Persister(path)

def _makePath(filename=SERIALIZE_FILE):
  this_dir = os.path.dirname(os.path.abspath(__file__))
  return os.path.join(this_dir, filename)

def _getData():
  trinary = TrinaryData(is_averaged=False,
      is_dropT1=False)
  return trinary.df_X, trinary.ser_y

def run(path, is_restart, max_iter=None,
    is_report=True):
  """
  Runs feature selection.
  :param str path: path to pcl file
  :param bool is_restart: Start a new analysis
  """ 
  df_X, ser_y = _getData()
  if max_iter is None:
    max_iter = len(df_X.columns) 
  persister = _getPersister(path)
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
        bcfo_kwargs=bcfo_kwargs)
  #
  optimizer.fit(df_X, ser_y)

def report(path=None):
  """
  Reports
  :param bool is_restart: Start a new analysis
  """
  def prt(header, dct):
    edct = ExtendedDict(dct)
    print("\n%s\n" % header)
    print(str(edct))
  #
  persister = _getPersister(path)
  if persister.isExist():
    optimizer = persister.get()
    prt("\n**Features by state:\n",
        optimizer.feature_dct)
    prt("\n**Scores by state:\n",
        optimizer.score_dct)
    prt("\n**Best scores by state:\n",
        optimizer.best_score_dct)
  else:
    print("***Serialization file not found: %s" % path)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      description="Run Feature Optimizer")
  parser.add_argument("-f", type=str,
      help="Pickle file with state",
      default=SERIALIZE_FILE)
  parser.add_argument("--restart",
      action="store_true",
      help="Re-start the run from beginning",
      )
  parser.add_argument("--only_report",
      action="store_true",
      help="Just provide a report"
      )
  args = parser.parse_args()
  path = _makePath(args.f)
  if args.only_report:
    report(path=path)
  else:
    run(path=path, is_restart=args.restart)
    report(path=path)
