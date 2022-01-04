"""Constructs trinary feature vectors for sample data."""

"""
feature vector is a dataframe
    index: instance
    column: gene
    value: trinary
features - names of genes

averaging is done in log2 units.
"""

import common.constants as cn
from common import trinary_data
from common import util
from common import transform_data

import collections
import os
import pandas as pd

# Reference types
REF_TYPE_BIOREACTOR = "ref_type_bioreactor"  # Use bioreactor data as reference
REF_TYPE_SELF = "ref_type_self"  # Internal reference
REF_TYPE_POOLED = "ref_type_pooled"  # Pool the data to get a reference
COL_TIME = "time"
COL_REP = "rep"
COL_CONDITION = "condition"
###
# Define characteristics of the sample files
###
# File names
FILE_GALAGAN = "galagan_raw_hypoxia_ts.csv"
FILE_AM_MDM = "AM_MDM_Mtb_transcripts_DEseq.csv"
FILE_AW = "AW_plus_v_AW_neg_Mtb_transcripts_DEseq.csv"
FILE_RUSTAD = "rustad_hypoxia_dataset_GSE9331.csv"
FILE_GSE167232 = "GSE167232_mtb_transcriptome_counts_normalized_filtered.csv"
# Functions for selecting instances that are reference values for gene expression
REFSEL_FUNC_AM_MDM = lambda i: ("AM" in i) and (not "1" in i)
REFSEL_FUNC_AW = lambda i: ("neg" in i) and (not "1" in i)
REFSEL_FUNC_GALAGAN = lambda i: ("d1." in i) and ("rep1" not in i)
REFSEL_FUNC_GSE167232 = None
REFSEL_FUNC_RUSTAD = lambda i: ("_4hr_" in i) 
# TODO: remove this logic since using all replications to calculate reference
#   and any([r in i for r in ["rep1", "rep2", "rep3"]])
# Strings that identify conditions (which is a grouping of replicas).
# For progression data, the condition strings are times. They are
# ordered in sequence from early to late.
CONDITION_STRINGS_AM_MDM = ["AM", "MDM"]
CONDITION_STRINGS_AW = ["AW_plus", "AW_neg"]
CONDITION_STRINGS_GALAGAN = ["t0h", "d1", "d2", "d3", "d5", "d7", "d8"]
CONDITION_STRINGS_GSE167232 = ["TB_HIGH", "TB_LOW", "TB_AM", "TB_IM"]
CONDITION_STRINGS_RUSTAD = ["H37Rv_hypoxia_%s" % s for s
    in [ "4hr", "8hr", "12hr", "1day", "4day", "7day"]]
# Strings that name replicas
REPLICA_NAMES_AM_MDM = ["_1", "_3", "_4", "_5"]
REPLICA_NAMES_AW = ["_1", "_3", "_4"]
REPLICA_NAMES_GALAGAN = ["rep1", "rep2", "rep3"]
REPLICA_NAMES_GSE167232 = ["1", "2", "3"]
REPLICA_NAMES_RUSTAD = ["rep1", "rep2", "rep3"]
# Describes data in the file
# csv: CSV file
# log2: True if in log2 units
# nrml: Normalized for gene length and library size
# sel: Reference selection function used in REF_TYPE_SELF
# cnm: Strings that identify conditions or time
# rnm: Strings that name the replicas
SampleDescriptor = collections.namedtuple("SampleDescriptor",
    "csv log2 nrml sel cnm rnm")
SAMPLE_DESCRIPTOR_DCT = {
    "AM_MDM": SampleDescriptor(csv=FILE_AM_MDM, log2=False, nrml=True,
              sel=REFSEL_FUNC_AM_MDM,
              cnm=CONDITION_STRINGS_AM_MDM,
              rnm=REPLICA_NAMES_AM_MDM),
    "AW": SampleDescriptor(csv=FILE_AW, log2=True, nrml=True,
          sel=REFSEL_FUNC_AW,
          cnm=CONDITION_STRINGS_AW,
          rnm=REPLICA_NAMES_AW),
    "galagan": SampleDescriptor(csv=FILE_GALAGAN, log2=True, nrml=True,
               sel=REFSEL_FUNC_GALAGAN,
               cnm=CONDITION_STRINGS_GALAGAN,
               rnm=REPLICA_NAMES_GALAGAN),
    "rustad": SampleDescriptor(csv=FILE_RUSTAD, log2=True, nrml=True,
              sel=REFSEL_FUNC_RUSTAD,
              cnm=CONDITION_STRINGS_RUSTAD,
              rnm=REPLICA_NAMES_RUSTAD),
    "GSE167232": SampleDescriptor(csv=FILE_GSE167232, log2=False, nrml=True,
                 sel=REFSEL_FUNC_GSE167232,
                 cnm=CONDITION_STRINGS_GSE167232,
                 rnm=REPLICA_NAMES_GSE167232),
    }
SAMPLES = list(SAMPLE_DESCRIPTOR_DCT.keys())


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
      is_average=False,
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
    is_average: bool
        Replicas are averaged
    """
    self.is_regulator = is_regulator
    self.is_display_errors = is_display_errors
    self.ref_type = ref_type
    self.is_average = is_average
    # Feature vectors for the samples
    self.df_AM_MDM = None
    self.df_AW = None
    self.df_galagan = None
    self.df_rustad = None
    self.df_GSE167232 = None

  def getDataframeAttributeName(self, sample_name):
    """
    Provides the attribute name for the sample data frame.

    Parameters
    ----------
    sample_name: str

    Returns
    -------
    str
    """
    return "df_%s" % sample_name

  def getDataframe(self, sample_name):
    """
    Provides the dataframe for the sample.

    Parameters
    ----------
    sample_name: str

    Returns
    -------
    pd.DataFrame
    """
    attribute_name = self.getDataframeAttributeName(sample_name)
    return self.__getattribute__(attribute_name)

  @property
  def AM_MDM(self):
    raise RuntimeError("Use `df_AM_MDM`")

  @property
  def AW(self):
    raise RuntimeError("Use `df_AW`")

  @property
  def sherman(self):
    raise RuntimeError("Unsupported data`")

  @property
  def galagan(self):
    raise RuntimeError("Use `df_galagan`")

  @property
  def rustad(self):
    raise RuntimeError("Use `df_rustad`")

  @property
  def GSE167232(self):
    raise RuntimeError("Use `df_GSE167232`")

  @staticmethod
  def _calcRef(df, selFunc=None):
    return SampleData._calcRefFromIndices(df, selFunc)

  def keys(self):
    return SAMPLES

  def values(self):
    return [self[s] for s in SAMPLES]

  def __getitem__(self, sample_name):
    """
    Provides the DataFrame for a sample using the syntax
        data[<sample_name>]

    Parameters
    ----------
    sample_name: str

    Returns
    -------
    DataFrame
    -------
    """
    return self.getDataframe(sample_name)

  def _makeSortKey(self, sample_name, instance_name):
    """
    Creates a sort key based on the order of the replica names.

    Parameters
    ----------
    sample_name: str
    instance_name: str
        instance as named in the index
    
    Returns
    -------
    int
    """
    descriptor = SAMPLE_DESCRIPTOR_DCT[sample_name]
    condition_strings = descriptor.cnm
    results = [i for i, s in enumerate(condition_strings) if s in instance_name]
    if len(results) == 1:
      result = results[0]
    elif len(results) != 0:
      result = len(condition_strings)
    else:
      raise RuntimeError(
         "Could not find replica string for sample %s, instance %s"
          % (sample_name, instance_name))
    return result

  def initialize(self):
    """
    Construct the feature vectors for the samples.
    """
    # Iterate across all samples
    for sample_name, descriptor in SAMPLE_DESCRIPTOR_DCT.items():
      attribute_name = self.getDataframeAttributeName(sample_name)
      ###
      # Construct a data frame that is normalized for gene and library
      # and has log2 units
      ###
      # Select indices that are for the conditions/times considered
      df = transform_data.readGeneCSV(descriptor.csv).T
      sel = [ any([d in i for d in descriptor.cnm]) for i in df.index]
      df = df[sel]
      # Sort the instances and complete initial processing
      indices = sorted(df.index, key=lambda v: self._makeSortKey(sample_name, v))
      df.index = indices
      if not descriptor.nrml:
        raise RuntimeError("Do gene normalization for sample %s" % sample_name)
      if not descriptor.log2:
        df = util.convertToLog2(df)
      ###
      # Convert to trinary values. This takes into account the reference values
      # for gene expression
      ###
      if self.ref_type == REF_TYPE_BIOREACTOR:
        ser_ref = transform_data.makeBioreactorT0ReferenceData()
      elif self.ref_type == REF_TYPE_POOLED:
        ser_ref = df.mean(axis=0)
      elif self.ref_type == REF_TYPE_SELF:
        if descriptor.sel is None:
          print("***%s: no selection for reference type 'self'. Result is None."
              % sample_name)
          df = None
        else:
          ser_ref = self._calcRefFromIndices(df, descriptor.sel)
      else:
        raise RuntimeError("%s is an invalid reference type" % self.ref_type)
      if df is not None:
        ###
        # Average replicas if requested
        ###
        if self.is_average:
          df = self.averageReplicas(df, descriptor.cnm)
        ###
        # Convert to trinary values
        ###
        df = transform_data.calcTrinaryComparison(df.T, ser_ref,
            is_convert_log2=False).T
        ###
        # Restrict to regulators?
        ###
        if self.is_regulator:
          trinary_data.subsetToRegulators(df)
      #
      self.__setattr__(attribute_name, df)

  @staticmethod
  def averageReplicas(df, condition_names):
    """
    Constructs a dataframe in which columns are averaged for indices that
    are for the same condition.

    Parameters
    ----------
    df: pd.DataFrame
    condition_names: list-str

    Returns
    -------
    pd.DataFrame
    """
    sers = []
    for name in condition_names:
      indices = [i for i in df.index if name in i]
      sers.append(df.loc[indices, :].mean())
    df_result = pd.concat(sers, axis=1)
    df_result = df_result.T
    df_result.index = condition_names
    return df_result

  @staticmethod
  def _calcRefFromIndices(df, selrefFunc):
    """
    Calculates the reference values for genes by using a selected set of indices.

    Parameters
    ----------
    df: DataFrame
    selrefFunc: Function
        parameters: indice
        returns: bool

    Returns
    -------
    Series
    """
    ref_idxs = [i for i in df.index if selrefFunc(i)]
    df_ref = df.loc[ref_idxs, :]
    return df_ref.mean()

  def serialize(self, directory=cn.TRINARY_SAMPLES_DIR):
    """
    Creates data in trinary feature matrix.
    :param SampleData sample_data:
    """
    for source in SAMPLES:
      path = os.path.join(directory, "%s.csv" % source)
      trinary_data.serializeFeatureMatrix(self.getDataframe(source), path)
