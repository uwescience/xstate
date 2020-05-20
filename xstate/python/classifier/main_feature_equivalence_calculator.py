'''Runs FeatureEquivalenceCalculator for one class.'''


import common.constants as cn
from common.data_provider import DataProvider
from common.trinary_data import TrinaryData
from common_python.classifier import  \
    multi_classifier_feature_optimizer as mcfo
from common_python.classifier import  \
    feature_equivalence_calculator as feq
from common_python.util.persister import Persister

import argparse
import os
import pandas as pd


NUM_STATES = 6
DIR = os.path.dirname(os.path.abspath(__file__))
PERSISTER_PATH = os.path.join(DIR,
    "main_feature_equivalence_calculator.pcl")
PERSISTER_BASE = os.path.join(DIR,
    "main_feature_equivalence_calculator_%d.pcl")
PERSISTER_PATHS = [PERSISTER_BASE % s
    for s in range(NUM_STATES)]
FIT_RESULT_PATH = os.path.join(cn.DATA_DIR,
    "fit_result_tf.xlsx")
# CSV file written
OUT_PATH = os.path.join(DIR,
    "main_feature_equivalence_calculator.csv")
OUT_BASE = os.path.join(DIR,
    "main_feature_equivalence_calculator_%d.csv")
OUT_PATHS = [OUT_BASE % s for s in range(NUM_STATES)]
NUM_CROSS_ITER = 150  # Cross validation iterations
MIN_SCORE = 0.9
CSV = "csv"
XLSX = "xlsx"


def getPersister(path=PERSISTER_PATH):
  return Persister(path)

def _getData(state):
  """
  Obtains data for a binary classifier for the class.
  :param int state: state for which classification is done
  :param pd.DataFrame, pd.Series:
  """
  provider = DataProvider()
  provider.do()
  trinary = TrinaryData(is_averaged=False,
      is_dropT1=False)
  columns = set(trinary.df_X.columns).intersection(
      provider.tfs)
  columns = list(columns)
  ser_y = trinary.ser_y.apply(lambda v:
    1 if v == state else 0)
  return trinary.df_X[columns], ser_y

def _remove(path):
  if os.path.isfile(path):
    os.remove(path)

def makeFitResult(state,
    fit_result_path=FIT_RESULT_PATH,
    min_score=MIN_SCORE):
  """
  Creates FitResults for a state.
  :param int state:
  :param str fit_result_path: Fit results
  :return list-mcfo.FitResult:
  """
  ext = os.path.split(fit_result_path)
  ext =ext[1].split(".")
  if ext[0] == CSV:
    df = pd.read_csv(fit_result_path)
  elif ext[1] == XLSX:
    df = pd.read_excel(fit_result_path)
  else:
    raise ValueError("Unsupported file type.")
  # Process groups
  dfg_dct = df.groupby([cn.STATE, cn.GROUP]).groups
  fit_results = []
  for key, indices in dfg_dct.items():
    this_state, this_group = key
    if this_state != state:
      continue
    sels = df.loc[indices, cn.FEATURE].tolist()
    group = df.loc[indices[0], cn.GROUP]
    score = df.loc[indices[0], cn.SCORE]
    if score >= min_score:
      fit_result = mcfo.FitResult(idx=group,
          sels=sels, sels_score=score)
      fit_results.append(fit_result)
  #
  return fit_results

def run(state, persister_path=PERSISTER_PATH,
    out_path=OUT_PATH,
    fit_results=None,
    is_restart=True, is_report=True,
    min_score=MIN_SCORE,
    **kwargs):
  """
  Runs feature selection.
  :param int state: State being analyzed
  :param str persister_path: path to pcl file
  :param str out_file: file where csv results are written
  :param bool is_restart: Start a new analysis
  :param dict kwargs: optional arguments for
      FeatureEquivalenceCalculator
  """ 
  # Initializations
  df_X, ser_y = _getData(state)
  persister = getPersister(persister_path)
  if is_restart:
    _remove(persister_path)
    _remove(out_path)
  if fit_results is None:
    fit_results = makeFitResult(state,
        fit_result_path=FIT_RESULT_PATH,
        min_score=min_score)
  else:
    fit_results = [f for f in fit_results
        if f.sels_score >= min_score]
  # Recover a previous run, if it exiss
  if persister.isExist():
    if is_report:
      print ("\n***PCL found: %s\n" % persister_path)
    calculator = persister.get()
  else:
    if is_report:
      print ("\n***Creating new PCL file: %s\n"
          % persister_path)
    calculator = feq.FeatureEquivalenceCalculator(
        df_X, ser_y,
        is_restart=is_restart,
        persister=persister, **kwargs)
  # Run the calculations
  calculator.run(fit_results)
  df = calculator.ria_dct[state]
  df[cn.STATE] = state
  df.to_csv(out_path)


if __name__ == '__main__':
  msg = "Run FeatureEquivalenceCalculator for state."
  parser = argparse.ArgumentParser(description=msg)
  parser.add_argument("state", 
      help="Expression state to evaluate",
      type=int)
  parser.add_argument("--restart",
      action="store_true",
      help="Re-start the run from beginning",
      )
  args = parser.parse_args()
  run(args.state, 
      persister_path=PERSISTER_PATHS[args.state], 
      out_path=OUT_PATHS[args.state],
      is_restart=args.restart, is_report=True)
