'''Assesses instances using cases'''
# TODO: Add gene description
# TODO: Add output file

import common.constants as cn
import common_python.constants as ccn
from tools.shared_data import SharedData

import argparse
import collections
import multiprocessing
import numpy as np
import os
import pandas as pd


INSTANCE = "instance"
MAX_SL = 0.001
REPORT_INTERVAL = 25  # Computations between reports
NUM_FSET = 100  # Number of feature sets examined
Arguments = collections.namedtuple("Arguments",
    "state df num_fset")


def _runState(arguments):
  """
  Does case evaluation for all instances for a single state.
  Run in multiple proceses concurrently.

  Parameters
  ----------
  state: int
  df_instance: pd.DataFrame
      Instances of feature vectors
  num_fset: int

  Return
  ------
  pd.DataFrame
      Columns: FEATURE_VECTOR, SIGLVL, STATE, INSTANCE
  """
  state = arguments.state
  df_instance = arguments.df
  num_fset = arguments.num_fset
  #
  shared_data = SharedData()
  fset_selector = lambda f: True
  dfs = []
  for instance in df_instance.index:
      ser_X = df_instance.loc[instance, :]
      collection = shared_data.collection_dct[state]
      df = collection.getEvaluationDF(ser_X,
          fset_selector=fset_selector, num_fset=num_fset,
          max_sl=MAX_SL)
      if len(df) == 0:
        df = pd.DataFrame({
            ccn.FEATURE_VECTOR: [np.nan],
            ccn.SIGLVL: [np.nan],
            })
      df[cn.STATE] = state
      df[INSTANCE] = instance
      dfs.append(df)
  return pd.concat(dfs)


def run(csv_handle, out_filename="report.csv",
    num_fset=NUM_FSET):
    """
    Processes the

    Parameters
    ----------
    ser_X: pd.DataFrame
        Feature vector for a single instance
        
    Returns
    -------
    None.
    """
    # Initializations
    out_path = os.path.join(cn.SAMPLES_DIR, out_filename)
    df_instance = pd.read_csv(csv_handle)
    # Iterate across instances
    shared_data = SharedData()
    num_process = len(shared_data.states)
    arguments_list = [
        Arguments(state=s, df=df_instance,
        num_fset=num_fset)
        for s in shared_data.states]
    with multiprocessing.Pool(num_process) as pool:
      results = pool.map(_runState, arguments_list)
    pool.join()
    df_report = pd.concat(results)
    df_report.to_csv(out_path, index=False)


if __name__ == '__main__':
  msg = """
Run state evaluations.
csv_file should be structured as:
   Rows are genes.
   Columns are instances.
   First column names the instances.
   Values are trinary.
"""
  parser = argparse.ArgumentParser(description=msg)
  parser.add_argument("csv_file",
      help="Differential exxpression data",
      type=argparse.FileType('r', encoding='UTF-8'))
  args = parser.parse_args()
  run(args.csv_file)
