'''Runs MultiClassifier and reports results.'''
"""
MultiClassifier is run with the following hyperparameters
  num_iteration: number of genes
  feature_selector: default
  max_degrade: 0.01 (difference from using all features)
"""


from common_python.classifier.multi_classifier  \
    import MultiClassifier
from common_python.types.extended_dict  \
    import ExtendedDict
from common_python.classifier import feature_selector
from common_python.util.persister import Persister
from common.trinary_data import TrinaryData

import argparse
import os


SERIALIZE_FILE = "main_multi_classifier.pcl"
MAX_DEGRADE = 0.01


def _makePath(filename=SERIALIZE_FILE):
  this_dir = os.path.abspath(os.path.dirname("__file__"))
  return os.path.join(this_dir, filename)

def _getData():
  trinary = TrinaryData(is_averaged=False,
      is_dropT1=False)
  return trinary.df_X, trinary.ser_y

def run(path, is_start, max_iter=None):
  """
  Runs feature selection.
  :param str path: path to pcl file
  :param bool is_start: Start a new analysis
  """ 
  df_X, ser_y = _getData()
  if max_iter is None:
    max_iter = len(df_X.columns) 
  persister = Persister(path)
  if is_start:
    if os.path.isfile(path):
      os.remove(path)
  if persister.isExist():
    print ("\n***Previous state found: %s\n" % path)
    multi_clf = persister.get()
  else:
    print ("\n***Creating new state: %s\n"
          % path)
    multi_clf = MultiClassifier(feature_selector_cls=  \
        feature_selector.FeatureSelector,
        max_iter=max_iter,
        max_degrade=MAX_DEGRADE)
  #
  multi_clf.fit(df_X, ser_y, persister=persister)

def report(path=None):
  """
  Reports
  :param bool is_start: Start a new analysis
  """
  if path is None:
    path = _makePath()
  multi_clf = MultiClassifier.getClassifier(
      path=path)
  if multi_clf is not None:
    dct = ExtendedDict(multi_clf.selector.feature_dct)
    print("\n**Features by state:\n")
    print(str(dct))
  else:
    print("***Serialization file not found: %s" % path)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      description="Run MultiClassifier")
  parser.add_argument("--f", type=str,
      help="Pickle file with classifier state",
      default=SERIALIZE_FILE)
  parser.add_argument("-restart",
      action="store_true",
      help="Re-start the run from beginning",
      )
  parser.add_argument("-only_report",
      action="store_true",
      help="Just provide a report"
      )
  args = parser.parse_args()
  path = _makePath(args.f)
  if args.only_report:
    report(path=path)
  else:
    run(path=path, is_start=args.restart)
    report(path=path)
