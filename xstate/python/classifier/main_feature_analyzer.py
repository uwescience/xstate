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
PERSISTER_PATH_PAT = os.path.join(DIR,
    "persister_%d.pcl")
IS_RESTART = True


def _getData(state, columns=None, **kwargs):
  """
  Obtains data for a binary classifier for the class.
  :param int state: state for which classification is done
  :param kwargs: dict: Options for TrinaryData
  :returns pd.DataFrame, pd.Series:
  """
  trinary = TrinaryData(**kwargs)
  ser_y = trinary.ser_y.apply(lambda v:
    1 if v == state else 0)
  if columns is None:
    df_X = trinary.df_X
  else:
    df_X = trinary.df_X[columns].copy()
  return df_X, ser_y

def run(state, out_dir_pat=OUT_PATH_DIR_PAT, num_cross_iter=NUM_CROSS_ITER,
    report_interval=REPORT_INTERVAL,
    columns=None, is_restart=IS_RESTART, **kwargs):
  """
  Runs feature selection.
  :param int state: State being analyzed
  :param dict kwargs: optional arguments for
       FeatureAnalyzer
  :param list-str columns: columns of df_X to use
  :param kwargs dict: arguments for TrinaryData
  """
  # Initializations
  df_X, ser_y = _getData(state, columns, **kwargs)
  analyzer = feature_analyzer.FeatureAnalyzer(
      CLF, df_X, ser_y, max_features_for_pairing=100,
      persister_path=PERSISTER_PATH_PAT % state,
      num_cross_iter=num_cross_iter, report_interval=report_interval)
  out_dir = out_dir_pat  % state
  _ = analyzer.serialize(out_dir,
      is_restart=is_restart)


if __name__ == '__main__':
  msg = "Run FeatureAnalyzer for state."
  parser = argparse.ArgumentParser(description=msg)
  parser.add_argument("state",
      help="Expression state to evaluate",
      type=int)
  parser.add_argument("--restart",
      action="store_true",
      help="Re-start the run from beginning",
      default=True,
      )
  args = parser.parse_args()
  run(args.state, num_cross_iter=NUM_CROSS_ITER,
      is_restart=args.restart,
      report_interval=REPORT_INTERVAL,
      is_regulator=False, is_dropT1=False, is_averaged=False)
