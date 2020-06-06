'''Runs FeatureAnalyzer for one class.'''


from common.trinary_data import TrinaryData
from common_python.classifier import feature_analyzer
import argparse
import os
from sklearn import svm


NUM_STATES = 6
DIR = os.path.dirname(os.path.abspath(__file__))
# Metrics calculated
OUT_PATH_DIR_PAT = os.path.join(DIR,
    "feature_analyzer_%d")
NUM_CROSS_ITER = 300  # Cross validation iterations
CLF = svm.LinearSVC()
REPORT_INTERVAL = 25  # Computations between reports
ANALYZE_METRICS = [feature_analyzer.SFA,
    feature_analyzer.CPC, feature_analyzer.IPA]
ANALYZE_METRICS = [feature_analyzer.IPA]


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

def run(state, out_dir_pat=OUT_PATH_DIR_PAT,
    columns=None, **kwargs):
  """
  Runs feature selection.
  :param int state: State being analyzed
  :param dict kwargs: optional arguments for
       FeatureAnalyzer
  :param list-str columns: columns of df_X to use
  """
  # Initializations
  df_X, ser_y = _getData(state, columns)
  analyzer = feature_analyzer.FeatureAnalyzer(
      CLF, df_X, ser_y, **kwargs)
  out_dir = out_dir_pat  % state
  analyzer.serialize(out_dir)


if __name__ == '__main__':
  msg = "Run FeatureAnalyzer for state."
  parser = argparse.ArgumentParser(description=msg)
  parser.add_argument("state",
      help="Expression state to evaluate",
      type=int)
  args = parser.parse_args()
  run(args.state, num_cross_iter=NUM_CROSS_ITER,
      report_interval=REPORT_INTERVAL)
