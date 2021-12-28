"""Constructs trinary feature vectors for sample data."""

"""
feature vector is a dataframe
    index: instance
    column: gene
    value: trinary
features - names of genes
"""

import common.constants as cn
from common import trinary_data
from common import util
from common.msg import writeMessage
from common.data_provider import DataProvider
from common import data_provider
import common.transform_data as transform_data

import collections
import copy
import numpy as np
import os
import pandas as pd

T1_INDEX = "T1"
MIN_NUM_NORMOXIA = 2  # Minimum number of normoxia states
FILE_GALAGAN = "galagan_raw_hypoxia_ts.csv" 
FILE_AM_MDM = "AM_MDM_Mtb_transcripts_DEseq.csv"
FILE_AW = "AW_plus_v_AW_neg_Mtb_transcripts_DEseq.csv"
FILE_RUSTAD = "rustad_hypoxia_dataset_GSE9331.csv"
FILE_GSE167232 = "GSE167232_mtb_transcriptome_counts_normalized_filtered.csv"
SHERMAN_INDUCED_PATH = os.path.join(cn.SAMPLES_DIR,
    "sherman_induced_mtb.txt")
SHERMAN_REPRESSED_PATH = os.path.join(cn.SAMPLES_DIR,
    "sherman_repressed_mtb.txt")
PROVIDER = DataProvider(is_display_errors=False)
PROVIDER.do()
# Reference types
REF_TYPE_BIOREACTOR = "ref_type_bioreactor"  # Use bioreactor data as reference
REF_TYPE_SELF = "ref_type_self"  # Internal reference
REF_TYPE_POOLED = "ref_type_pooled"  # Pool the data to get a reference
COL_TIME = "time"
COL_REP = "rep"
COL_CONDITION = "condition"
SAMPLES = ["AM_MDM", "AW", "sherman", "galagan", "rustad", "GSE167232"]


################## FUNCTIONS ###############
def getSampleData(**kwargs):
  data = SampleData(**kwargs)
  data.initialize()
  return data


################## CLASSES ###############
class SampleData(object):

  def __init__(self, is_regulator=True,
      is_display_errors=False,
      ref_type=REF_TYPE_BIOREACTOR,
    ):
    """
    Acquires data obtain from other soruces.
    
    Parameters
    ----------
    is_regulator: bool
        Only return gene regulators
    is_display_error: bool
        Report errors in constructing sample data
    ref_type: str
        What reference data are used to calculate gene expression
    """
    self.is_regulator = is_regulator
    self.is_display_errors = is_display_errors
    self.ref_type = ref_type
    # Feature vectors for the samples
    self.df_AM_MDM = None
    self.df_AW = None
    self.df_sherman = None
    self.df_galagan = None
    self.df_rustad = None
    self.df_GSE167232 = None

  def getDataFrame(self, sample_name):
    """
    Provides the dataframe for the sample.

    Parameters
    ----------
    sample_name: str
    
    Returns
    -------
    pd.DataFrame
    """
    attribute_name = "df_%s" % sample_name
    return self.__getattribute__(attribute_name)

  @property
  def AM_MDM(self):
    raise RuntimeError("Use `df_AM_MDM`")

  @property
  def AW(self):
    raise RuntimeError("Use `df_AW`")

  @property
  def sherman(self):
    raise RuntimeError("Use `df_sherman`")

  @property
  def galagan(self):
    raise RuntimeError("Use `df_galagan`")

  @property
  def rustad(self):
    raise RuntimeError("Use `df_rustad`")

  @property
  def GSE167232(self):
    raise RuntimeError("Use `df_GSE167232`")

  def initialize(self):
    """
    Construct the feature vectors for the samples.
    """
    ####
    # AM/MDM
    ####
    if (self.ref_type == REF_TYPE_BIOREACTOR):
      self.df_AM_MDM = transform_data.trinaryReadsDF(
          is_display_errors=self.is_display_errors,
          csv_file=FILE_AM_MDM,
          is_time_columns=False).T
    elif (self.ref_type == REF_TYPE_SELF):
      ref_sel_func = lambda i: ("AM" in i) and (not "1" in i)
      def calcRef(df):
        return self._calcRefFromIndices(df, ref_sel_func)
      self.df_AM_MDM = self._makeSamplesWithPooledReference(FILE_AM_MDM, calcRef, False)
    else:
      self.df_AM_MDM = self._makeSamplesWithPooledReference(FILE_AM_MDM, self.calcRefPooled, False)
    ####
    # AW
    ####
    if (self.ref_type == REF_TYPE_BIOREACTOR):
      self.df_AW = transform_data.trinaryReadsDF(
          csv_file=FILE_AW,
          is_display_errors=self.is_display_errors,
          is_time_columns=False).T
    elif (self.ref_type == REF_TYPE_SELF):
      ref_sel_func = lambda i: ("neg" in i) and (not "1" in i)
      def calcRef(df):
        return self._calcRefFromIndices(df, ref_sel_func)
      self.df_AW = self._makeSamplesWithPooledReference(FILE_AW, calcRef, False)
    else: # Pooled
      self.df_AW = self._makeSamplesWithPooledReference(FILE_AW, self.calcRefPooled, False)
    self.df_AW = self.df_AW.sort_index()
    ####
    # Galagn data
    ####
    if (self.ref_type == REF_TYPE_BIOREACTOR):
      self.df_galagan = self._getGalaganData()
      self.df_galagan = trinary_data.convertToTrinary(self.df_galagan)
    elif (self.ref_type == REF_TYPE_SELF):
      df_data = self._getGalaganData()
      ref_sel_func = lambda i: ("d1." in i) and ("rep1" not in i)
      def calcRef(df):
        return self._calcRefFromIndices(df, ref_sel_func)
      self.df_galagan = self._makeSamplesWithPooledReference(
          None, calcRef, False, df_data=df_data)
    else:  # Pooled
      df_data = self._getGalaganData()
      self.df_galagan = self._makeSamplesWithPooledReference(
          None, self.calcRefPooled, False, df_data=df_data)
    ####
    # Sherman
    ####
    self.df_sherman = trinary_data.getTrinaryFromGeneLists()
    self.df_sherman = self.df_sherman.transpose()
    ####
    # Rustad
    ####
    if (self.ref_type == REF_TYPE_BIOREACTOR):
      self.df_rustad = transform_data.trinaryReadsDF(
          csv_file=FILE_RUSTAD,
          is_display_errors=self.is_display_errors,
          is_convert_log2=False,
          is_time_columns=False).T
    elif (self.ref_type == REF_TYPE_SELF):
      ref_sel_func = lambda i: ("_4hr_" in i) and ("rep6" not in i)
      def calcRef(df):
        return self._calcRefFromIndices(df, ref_sel_func)
      self.df_rustad = self._makeSamplesWithPooledReference(FILE_RUSTAD, calcRef, False)
    else: # pooled
      self.df_rustad = self._makeSamplesWithPooledReference(FILE_RUSTAD, self.calcRefPooled, False)
    # Construct the major sor index for rustad
    time_vals = ["4hr", "8hr", "12hr", "1day", "4day", "7day"]
    reps  = [i.split("_")[-1] for i in self.df_rustad.index]
    times  = [i.split("_")[-2] for i in self.df_rustad.index]
    conditions = ["_".join(i.split("_")[0:2]) for i in self.df_rustad.index]
    self.df_rustad[COL_TIME] = [time_vals.index(v) for v in times]
    self.df_rustad[COL_REP] = reps
    self.df_rustad[COL_CONDITION] = conditions
    self.df_rustad = self.df_rustad.sort_values([COL_CONDITION, COL_TIME, COL_REP])
    for col in [COL_REP, COL_TIME, COL_CONDITION]:
      del self.df_rustad[col]
    ####
    # GSE167232
    ####
    if (self.ref_type == REF_TYPE_BIOREACTOR)  \
        or (self.ref_type == REF_TYPE_SELF):
      if (self.ref_type == REF_TYPE_SELF):
        message = "\n**No self reference defined for GSE167232."
        message += " Using bioreactor data.\n"
        writeMessage(message)
      self.df_GSE167232 = transform_data.trinaryReadsDF(
          csv_file=FILE_GSE167232,
          is_display_errors=self.is_display_errors,
          is_normalized=True,
          is_time_columns=False).T
    else: # pooled
      self.df_GSE167232 = self._makeSamplesWithPooledReference(FILE_GSE167232,
          self.calcRefPooled, False)
    ###
    # Restrict to regulators?
    ###
    if self.is_regulator:
      for df in [self.df_AM_MDM, self.df_AW, self.df_sherman, self.df_galagan,
          self.df_rustad, self.df_GSE167232]:
        trinary_data.subsetToRegulators(df)

  @staticmethod
  def _makeSamplesWithPooledReference(csv_file, calcRef, is_time_columns,
      df_data=None):
    """
    Constructs samples with respect to reference instances within the same
    data set.

    Parameters
    ----------
    csv_file: str
        CSV file for the sample data if df_data == None
    calcRef: Function
        parameters: DataFrame
        returns: Series
    is_time_columns: bool
        has a "time" column
    df_data: dataframe
        count (not log2) to be converted into TrinaryData
    Returns
    -------
    DataFrame: Trinary values of samples
    """
    if df_data is None:
      df_data = transform_data.readGeneCSV(csv_dir=cn.SAMPLES_DIR,
          csv_file=csv_file).T
    df_normalized = PROVIDER.normalizeReadsDF(
        df_data.T, is_time_columns=False).T
    df_normalized = df_normalized.applymap(lambda v: max(v, cn.MIN_VALUE))
    ser_ref = calcRef(df_normalized)
    ser_ref = ser_ref.map(lambda v: max(v, cn.MIN_VALUE))  # No 0's
    ser_ref_log = util.convertToLog2(ser_ref)
    df_normalized_log = util.convertToLog2(df_normalized)
    df_sample_trinary = transform_data.calcTrinaryComparison(
        df_normalized_log.T,
        ser_ref=ser_ref_log, is_convert_log2=False).T
    return df_sample_trinary

  @staticmethod
  def _calcRefFromIndices(df, sel_ref_func):
    ref_idxs = [i for i in df.index if sel_ref_func(i)]
    df_ref = df.loc[ref_idxs, :]
    return df_ref.mean()

  @staticmethod
  def calcRefPooled(df):
    return df.mean(axis=0)

  def _getGSE167232(self):
    df = transform_data.trinaryReadsDF(
        csv_file=FILE_GSE167232,
        is_display_errors=self.is_display_errors,
        is_normalized=True,
        is_time_columns=False).T
    if self.is_regulator:
      trinary_data.subsetToRegulators(df)
    return df
  
  def _getGalaganData(self):
    """
    Constructs trinary values for Galagan data.
    These data are normalized and in log2 units.
    The 10h data are Normoxia.
    :return pd.DataFrame:
        columns: genes
        index: time instances
        values: trinary based on log2
    """
    self.df_galagan = transform_data.readGeneCSV(
        csv_file=FILE_GALAGAN)
    dfs = []
    for idx in range(1, 4):
      stg = "rep%d" % idx
      columns = [c for c in self.df_galagan.columns
          if stg in c]
      col_ref = columns[0]
      columns.remove(col_ref)
      df = self.df_galagan[columns].copy()
      df = df.apply(lambda c: c - self.df_galagan[col_ref])
      dfs.append(df)
    df_merge = pd.concat(dfs, axis=1)
    return df_merge.T
