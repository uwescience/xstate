"""Classifies expression data."""

import add_path
from common import constants as cn
from common import trinary_data
from common import transform_data
from tools import make_svm_ensemble

import argparse
import os
import numpy as np
import pandas as pd


OUTPUT_PATH = "classify_output.csv"


def processExpressionData(data_file, data_dir=cn.SAMPLES_DIR,
    parameter_file=cn.ENSEMBLE_PATH,
    output_path=OUTPUT_PATH, is_display_errors=True):
  """
  Performs classification on expression data.
  Expression data is structured as a comma separated variable (CSV)
  file. Rows are instances to classify. Columns are genes.
  The first row (header row) contains the names of genes for columns.
  :param str data_dir: directory containing the data samples
  :param str data_path: path to the expression data
  :param str parameter_file: path containing ensemble parameters
  Writes classifications to standard out.
  """
  svm_ensemble = make_svm_ensemble.make(file_path=parameter_file,
      is_force=False)
  trinary = trinary_data.TrinaryData(
      is_display_errors=is_display_errors)
  inv_dict = {v: k for k, v in trinary.state_dict.items()}
  path = os.path.join(data_dir, data_file)
  df_trinary = transform_data.trinaryReadsDF(csv_file=path,
      is_display_errors=is_display_errors)
  df_classes = svm_ensemble.predict(df_trinary.T)
  columns = [inv_dict[v] for v in df_classes.columns]
  df_classes.columns = [inv_dict[v] for v in df_classes.columns]
  df_classes.to_csv(output_path)


if __name__ == '__main__':
  # Do arg parse with errors
  desc = 'Classifies gene expression data.'
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('dfile',
      help='file containing expression data (csv)',
      type=str)
  parser.add_argument('--dir', '-d',
       help='directory containing input expression data files',
      default=cn.SAMPLES_DIR)
  parser.add_argument('--pfile', '-p',
       help='classifier parameter file',
      default=cn.ENSEMBLE_PATH)
  parser.add_argument('--ofile', '-o',
       help='output file path with classifier results (csv)',
      default=OUTPUT_PATH)
  args = parser.parse_args()
  processExpressionData(args.dfile, 
      data_dir=args.dir, parameter_file=args.pfile,
      output_path=args.ofile)
  print("Success!")
