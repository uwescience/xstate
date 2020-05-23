'''Runs FeatureAnalyzer for one class.'''


import common.constants as cn
from common.trinary_data import TrinaryData
from common_python.classifier import feature_analyzer
import argparse
import os
import pandas as pd
from sklearn import svm


NUM_STATES = 6
DIR = os.path.dirname(os.path.abspath(__file__))
# Metrics calculated
OUT_PATH_PAT = os.path.join(DIR,
    "main_feature_analyzer_%s_%d.csv")
NUM_CROSS_ITER = 150  # Cross validation iterations
CLF = svm.LinearSVC()
REPORT_INTERVAL = 25  # Computations between reports


def _getData(state, columns=None):
  """
  Obtains data for a binary classifier for the class.
  :param int state: state for which classification is done
  :param pd.DataFrame, pd.Series:
  """
  trinary = TrinaryData(is_averaged=False,
      is_dropT1=False, is_regulator=True)
  ser_y = trinary.ser_y.apply(lambda v:
    1 if v == state else 0)
  if columns is None:
    df_X = trinary.df_X
  else:
    df_X = trinary.df_X[columns].copy()
  return df_X, ser_y

def _msg(metric, state, is_report):
  if is_report:
    print ("\n***Completed metric %s for state %d\n" %
        (metric, state))

def run(state, out_path_pat=OUT_PATH_PAT, is_report=True,
    columns=None, **kwargs):
  """
  Runs feature selection.
  :param int state: State being analyzed
  :param Pattern out_path_path: pattern for output files
      _%s: metric
      _%d: state
  :param dict kwargs: optional arguments for
       FeatureAnalyzer
  :param list-str columns: columns of df_X to use
  """ 
  # Initializations
  df_X, ser_y = _getData(state, columns)
  analyzer = feature_analyzer.FeatureAnalyzer(
      CLF, df_X, ser_y, **kwargs)
  #
  def analyze(metric):
    path = out_path_pat % (metric, state)
    if metric == feature_analyzer.SFA:
      result = analyzer.ser_sfa
    elif metric == feature_analyzer.CPC:
      result = analyzer.df_cpc
    elif metric == feature_analyzer.IPA:
      result = analyzer.df_ipa
    result.to_csv(path)
    _msg(metric, state, is_report)
  # Process
  for metric in feature_analyzer.METRICS:
    analyze(metric)


if __name__ == '__main__':
  msg = "Run FeatureAnalyzer for state."
  parser = argparse.ArgumentParser(description=msg)
  parser.add_argument("state", 
      help="Expression state to evaluate",
      type=int)
  args = parser.parse_args()
  run(args.state, num_cross_iter=NUM_CROSS_ITER,
      report_interval=REPORT_INTERVAL)
