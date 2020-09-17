'''Assesses instances using cases'''
# TODO: Add output file

import common.constants as cn
from common.data_provider import DataProvider
import common_python.constants as ccn
from common_python.classifier.feature_set import  \
    FeatureVector
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
INDEX = "index"


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
      FEATURE_VECTOR
      SIGLVL: significance level of FRAC
      STATE: state analyzed
      INSTANCE: from data feature vector
      COUNT: number of cases
      FRAC: fraction of positive cases
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
      df = collection.getFVEvaluations(ser_X,
          fset_selector=fset_selector, num_fset=num_fset,
          max_sl=MAX_SL)
      if len(df) > 0:
        df[cn.STATE] = state
        df[INSTANCE] = instance
      dfs.append(df)
  df_result = pd.concat(dfs)
  df_result.index = range(len(df_result.index))
  # Augment the dataframe with gene descriptions
  provider = DataProvider()
  provider.do()
  df_go = provider.df_go_terms
  descriptions = []
  for stg in df_result[ccn.FEATURE_VECTOR]:
    if not isinstance(stg, str):
      descriptions.append("")
    else:
      feature_vector = FeatureVector.make(stg)
      features = feature_vector.fset.set
      description = []
      for feature in features:
        df_sub = df_go[df_go[cn.GENE_ID] == feature]
        this_desc = ["%s: %s " % (feature, f)
            for f in df_sub[cn.GO_TERM]]
        description.extend(this_desc)
      description = "\n".join(description)
      descriptions.append(description)
  #
  df_result[cn.GENE_DESCRIPTION] = descriptions
  return df_result


def run(input_fd, output_fd, num_fset=NUM_FSET):
    """
    Processes the

    Parameters
    ----------
    input_fd: File Descriptor
        Input CSV file
    output_fd: File Descriptor
        Output file
    num_fset: int
        Number of FeatureSets considered
        
    Returns
    -------
    None.
    """
    # Initializations
    df_instance = pd.read_csv(input_fd)
    if not INDEX in df_instance.columns:
       msg = "One input column must be named 'instance'."
       raise ValueError(msg)
    df_instance = df_instance.set_index(INDEX)
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
    df_report.to_csv(output_fd, index=False)
    input_fd.close()
    output_fd.close()


if __name__ == '__main__':
  msg = """
Run state evaluations.
csv_file should be structured as:
   Rows are genes.
   Columns are instances. There should be an 'index' column name
   First column names the instances.
   Values are trinary.
"""
  parser = argparse.ArgumentParser(description=msg)
  parser.add_argument("input_file",
      help="Input CSV file of trinary exxpression data",
      type=argparse.FileType('r', encoding='UTF-8'))
  parser.add_argument("output_file",
      help="Output CSV file of report by feature vector",
      type=argparse.FileType('w', encoding='UTF-8'))
  args = parser.parse_args()
  run(args.input_file, args.output_file)
