'''Runs FeatureAnalyzer for one class.'''


from common.trinary_data import TrinaryData
from common_python.classifier import feature_analyzer
from common_python.util.persister import Persister
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
MAX_FEATURES_FOR_PAIRING = 100


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
    is_status=False,
    report_interval=REPORT_INTERVAL,
    is_report=True,
    columns=None, is_restart=IS_RESTART, **kwargs):
  """
  Runs feature selection.
  :param int state: State being analyzed
  :param list-str columns: columns of df_X to use
  :param bool is_status: report status extracted from the persister
  :param kwargs dict: arguments for TrinaryData
  :param dict kwargs: optional arguments for
       FeatureAnalyzer
  """
  def calcLen(obj, func):
    if obj is None:
      return 0
    else:
      return func(obj)
  #
  CUR_LEN = "cur_len"
  MAX_LEN = "max_len"
  persister_path = PERSISTER_PATH_PAT % state
  df_X, ser_y = _getData(state, columns, **kwargs)
  if is_status:
    func = lambda d: len(d["score"])
    #
    persister = Persister(persister_path)
    analyzer = persister.get()
    pair_length = MAX_FEATURES_FOR_PAIRING*(MAX_FEATURES_FOR_PAIRING-1)/2  \
        + MAX_FEATURES_FOR_PAIRING
    dct = {
        "sfa": {CUR_LEN: calcLen(analyzer._sfa_dct, lambda d: len(d.keys())),
                MAX_LEN: len(df_X.columns)},
        "cpc": {CUR_LEN: calcLen(analyzer._cpc_dct, func),
                MAX_LEN: pair_length},
        "ipa": {CUR_LEN: calcLen(analyzer._ipa_dct, func),
                MAX_LEN: pair_length},
        }
    report_stg = "State %s: " % str(state)
    for metric in dct.keys():
      cur_length = 0
      if dct[metric][MAX_LEN] is not None:
        cur_length = dct[metric][CUR_LEN]
      frac = min(1.0, cur_length/dct[metric][MAX_LEN])
      report_stg =  ("%s %s/%2.3f" % (report_stg, metric, frac))
    if is_report:
      print(report_stg)
  else:
    analyzer = feature_analyzer.FeatureAnalyzer(
        CLF, df_X, ser_y,
        max_features_for_pairing=MAX_FEATURES_FOR_PAIRING,
        persister_path=persister_path,
        num_cross_iter=num_cross_iter, report_interval=report_interval)
    out_dir = out_dir_pat  % state
    _ = analyzer.serialize(out_dir, is_restart=is_restart)


if __name__ == '__main__':
  msg = "Run FeatureAnalyzer for state."
  parser = argparse.ArgumentParser(description=msg)
  parser.add_argument("state",
      help="Expression state to evaluate",
      type=int)
  parser.add_argument("--restart",
      action="store_true",
      help="Re-start the run from beginning",
      default=False,
      )
  parser.add_argument("--status",
      action="store_true",
      help="Report progress status for calculating feature metrics",
      default=False,
      )
  args = parser.parse_args()
  run(args.state, num_cross_iter=NUM_CROSS_ITER,
      is_restart=args.restart,
      is_status=args.status,
      report_interval=REPORT_INTERVAL,
      is_regulator=False, is_dropT1=False, is_averaged=False)
