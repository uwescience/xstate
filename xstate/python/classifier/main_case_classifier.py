'''Assesses instances using cases'''
# TODO: Multi-process the state analysis
# TODO: Profile

import common.constants as cn
from tools.shared_data import SharedData

import argparse
import numpy as np
import os
import pandas as pd


INSTANCE = "instance"
MAX_SL = 0.001
REPORT_INTERVAL = 25  # Computations between reports
DATA_PATH = os.path.join(cn.DATA_DIR, "feature_analyzer")


def run(csv_handle, out_filename="report.csv"):
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
    shared_data = SharedData()
    out_path = os.path.join(SAMPLE_DIR, out_filename)
    df_instance = pd.read_csv(csv_handle)
    fset_selector = lambda f: True
    # Iterate across instances
    dfs = []
    for idx, state in enumerate(shared_data.states):
      for instance in df_instance.index:
          collection = shared_data.collection_dct[state]
          df = collection.getEvaluationDF(ser_X,
              fset_selector=fset_selector, num_fset=100,
              MAX_SL=MAX_SL, **kwargs)
          df[cn.STATE] = state
          df[INSTANCE] = instance
          dfs.append(df)
    df_report = pd.concat(dfs)
    df_report.to_csv(out_path)


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
