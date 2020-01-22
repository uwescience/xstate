"""
Makes Data Available in Standard Formats. Creates the following data:
  df_gene_description DF
    cn.GENE_ID (index), cn.GENE_NAME, cn.LENGTH, cn.PRODUCT, cn.START, cn.END, cn.STRAND
  Dataframes in dfs_centered_adjusted_read_count
    cn.GENE_ID (index), columns: time indices
  hypoxia curve DF
    cn.SAMPLE, cn.HOURS, 0, 1, 2 (DO values), mean, std, cv
  df_normalized
     index: cn.GENE_ID
     column: time
  df_mean,   # Mean values of counts
  df_std,   # Std of count values
  df_cv,   # Coefficient of variation
  df_stage_matrix,
  df_gene_expression_state,   # Genes expressed in each state
  df_go_terms
    cn.GENE_ID cn.GO_TERM
  df_ec_terms
    cn.GENE_ID cpn.KEGG_EC
  df_ko_terms
    cn.GENE_ID cpn.KEGG_KO
  df_kegg_pathways
    cpn.KEGG_PATHWAY cpn.DESCRIPTION
  df_kegg_gene_pathways
    cn.GENE_ID cpn.KEGG_PATHWAY
  dfs_read_count - raw read count dataframes
     index: cn.GENE_ID
     column: time
  dfs_adjusted_read_count - readcounts adjusted w.r.t. library size, gene length
     index: cn.GENE_ID
     column: time
  dfs_adjusted_read_count_wrt0 - adjusted w.r.t. libary size, gene length, time 0
     index: cn.GENE_ID
     column: time
  dfs_centered_adjusted_read_count - centers w.r.t. mean value of gene
     index: cn.GENE_ID
     column: time
"""

import common.constants as cn
import common_python.constants as cpn
from common_python.util.persister import Persister

import os
import pandas as pd
import numpy as np


FILENAME_HYPOXIA = "hypoxia_curve_DO"
FILENAME_NORMALIZED = "normalized_log2_transformed_counts"
FILENAME_READS = "hypoxia_timecourse_reads"
FILENAME_STAGES = "stages_matrix"
FILENAME_GENEDATA = "gene_data"
FILENAME_GENE_EXPRESSION_STATE = "gene_expression_state"
FILENAME_GO_TERMS = "MTB.GO.All.GOterms"
FILENAME_EC_TERMS = "mtb_gene_ec"
FILENAME_KO_TERMS = "mtb_gene_ec"
FILENAME_KEGG_PATHWAYS = "mtb_kegg_pathways"
FILENAME_KEGG_GENE_PATHWAY = "mtb_kegg_gene_pathway"
NUM_REPL = 3
TIME_0 = "T0"
T0 = 0
MIN_LOG2_VALUE = -10
MIN_VALUE = 10e-5
MILLION = 1e6
KILOBASE = 1e3  # Thousand bases

class DataProvider(object):
  # Instance variables in the class
  instance_variables = [
    "df_gene_description",
    "df_gene_expression_state",
    "df_hypoxia",
    "df_mean",
    "df_std",
    "df_cv",
    "df_normalized",
    "df_stage_matrix",
    "df_go_terms",
    "df_ec_terms",
    "df_ko_terms",
    "df_kegg_pathways",
    "df_kegg_gene_pathways",
    "dfs_read_count",
    "_dfs_adjusted_read_count",
    "_dfs_adjusted_read_count_wrt0",
    "_dfs_centered_adjusted_read_count",
    ]

  def __init__(self, data_dir=cn.DATA_DIR, is_normalized_wrtT0=True,
      is_only_qgenes=True, is_display_errors=True):
    """
    :param bool is_normalized_wrtT0: normalize data w.r.t. T0
        Otherwise, standardize values using the mean.
    :param bool is_only_qgenes: only include genes included in multi-hypothesis test
    """
    self._data_dir = data_dir
    self._is_normalized_wrtT0 = is_normalized_wrtT0
    self._is_only_qgenes = is_only_qgenes
    self._is_display_errors = is_display_errors
    self._setValues()

  def _setValues(self, provider=None):
    """
    Sets values for the instance variables.
    :param DataProvider provider:
    """
    for var in self.__class__.instance_variables:
      if provider is None:
        stmt = "self.%s = None" % var
      else:
        stmt = "self.%s = provider.%s" % (var, var)
      exec(stmt)

  def _makeDFFromCSV(self, filename, is_index_geneid=False):
    """
    Processes a CSV file
    :param str filename: without csv extension
    :param bool is_index_geneid: use cn.GENE_ID to index
    :return pd.DataFrame:
    """
    path = os.path.join(self._data_dir, "%s.csv" % filename)
    df = pd.read_csv(path)
    if is_index_geneid:
      df = df.set_index(cn.GENE_ID)
    return df

  def _makeHypoxiaDF(self):
    """
    :return pd.DataFrame:
    Columns: HOURS, SAMPLE, MEAN, STD,
       0, 1, 2 are the dissolved oxygen for the replications
    """
    df = self._makeDFFromCSV(FILENAME_HYPOXIA)
    new_columns = {
        "sample": cn.SAMPLE,
        "time (h)":  cn.HOURS,
        "DO_reactor_A": 0,
        "DO_reactor_B": 1,
        "DO_reactor_C": 2,
        }
    df = df.rename(columns=new_columns)
    df_do = pd.DataFrame([df[n] for n in range(3)])
    df[cn.MEAN] = df_do.mean(axis=0)
    df[cn.STD] = df_do.std(axis=0)
    df[cn.CV] = 100 * df[cn.STD] / df[cn.MEAN]
    return df

  def _makeGeneDescriptionDF(self):
    df = self._makeDFFromCSV(FILENAME_GENEDATA)
    df[cn.LENGTH] = [np.abs(r[cn.START] - r[cn.END]) for _, r in df.iterrows()]
    df = df.rename(
        columns={"Locus Tag": cn.GENE_ID})
    return df.set_index(cn.GENE_ID)

  def _getNumRepl(self):
      return NUM_REPL

  def _makeMeanDF(self, is_abs=True):
      """
      Creates a dataframe for the mean values.
      :param bool is_abs: computes the mean of absolute values
      :return pd.DataFrame:
      """
      if is_abs:
        predicate = lambda v: np.abs(v)
      else:
        predicate = lambda v: v
      dfs = self.dfs_centered_adjusted_read_count
      dfs_new = [df.applymap(predicate) for df in dfs]
      return sum (dfs_new) / len(dfs)

  def _makeStdDF(self):
      """
      Creates a dataframe for the standard deviations.
      :return pd.DataFrame:
      """
      num_repl = self._getNumRepl()
      df_mean = self._makeMeanDF()
      dfs = self.dfs_centered_adjusted_read_count
      df_std = (sum([dfs[n]*dfs[n] for n in range(num_repl)])
          - num_repl * df_mean * df_mean) / (num_repl - 1)
      return df_std.pow(1./2)

  def _reduceDF(self, df):
    """
    Reduces a dataframe indexed by GENE_ID to those genes available in df_gene_description
    :param pd.DataFrame df: indexed by GENE_ID
    :return pd.DataFrame:
    """
    return pd.DataFrame([r for idx, r in df.iterrows() if idx in self.df_gene_description.index])

  def _makeReadCountDFS(self):
    """
    Creates a list of dataframes for each replication of the counts.
    :return list-pd.DataFrame:
      indexed by GENE_ID
      column names are integers of time indices
    Notes
      1. Assumes that self.df_gene_description has been constructed
      2. Counts are normalized:
         a. Adjusted for gene length
         b. Adjusted for library size of the replication
         c. Subtract the mean so that down regulation is negative
            and upregulation is positive.
    """
    dfs = []
    # Get the read counts and reduce it
    df_data = self._makeDFFromCSV(FILENAME_READS)
    df_data.index = df_data[cn.GENE_ID]
    df = self._reduceDF(df_data)
    # Separate the replications into different dataframes
    name_map = {0: 'A', 1: 'B', 2: 'C'}
    for repl in range(NUM_REPL):
        # Select the columns for this replication
        col_suffix = "_%s" % name_map[repl]
        column_names = [c for c in df.columns if col_suffix in c]
        # Transform names into numbers
        new_names = {}
        for name in column_names:
            split_name = name.split("_")
            new_name = split_name[2][1:]
            new_names[name] = int(new_name)
        df_repl = df[column_names]
        df_repl = df_repl.rename(columns=new_names)
        df_repl.index = df.index
        dfs.append(df_repl)
    #
    return dfs

  def _makeNormalizedDF(self):
    """
    Transformation of the "normalized Read Counts" processed by DESeq2.
    Standardized the values for each gene.
    Drops rows where all columns are minimum values.
    Assumes that self.df_gene_expression_state has been initialized.
    Only includes genes that are expressed.
    """
    df = self._makeDFFromCSV(FILENAME_NORMALIZED)
    df = df.set_index(cn.GENE_ID)
    # Normalize w.r.t. time 0
    if self._is_normalized_wrtT0:
      drops = []  # Rows to drop
      for idx in df.index:
        values = df.loc[idx, :] - df.loc[idx, TIME_0]
        df.loc[idx, :] = [max(MIN_LOG2_VALUE,  v) for v in values]
        if all([v <= MIN_LOG2_VALUE for v in df.loc[idx, :]]):
          drops.append(idx)
      df = df.drop(index=drops) # Drop the 0 rows
    # Normalize by standardizing the row
    else:
      for col in df.columns:
        df[col] = (df[col] - np.mean(df[col])) / np.std(df[col])
    # Find genes to keep
    if self._is_only_qgenes:
      keep_genes = self.df_gene_expression_state.index
      df = df[df.index.isin(keep_genes)]
    #
    return df
 
  def normalizeReadsDF(self, df):
    """
    Transforms a vector of read counts into features and units used in analysis.
      1. Adjusts read counts for each gene based on length and library size.
         Uses the metric RPKM - reads per kilobase millions
      2. Deletes unused features
    :param pd.DataFrame df: 
        index are genes, columns are instances, values are readcounts
    :return pd.DataFrame: 
        index are genes, columns are instances, values are readcounts
    """
    # Adjust for library size
    ser_tot = df.sum(axis=0)/MILLION
    df_result = df/ser_tot
    # Adjust for gene length
    for gene in df.index:
      if gene in self.df_gene_description.index:
        df_result.loc[gene, :] = KILOBASE*df_result.loc[gene,:] \
            / self.df_gene_description.loc[gene, cn.LENGTH]
      else:
        if self._is_display_errors:
          msg = "**Warning: Data doesn't have gene %s" % gene
          print(msg)
    # Find genes to keep
    if self._is_only_qgenes:
      keep_genes = self.df_gene_expression_state.index
      df_result = df_result[df_result.index.isin(keep_genes)]
    #
    return df_result

  def _makeStageMatrixDF(self):
    """
    Columns: STAGE_NAME, STAGE_COLOR
    Index: TIMEPOINT
    """
    df = self._makeDFFromCSV(FILENAME_STAGES)
    return df.set_index(cn.TIMEPOINT)

  def equals(self, provider):
    """
    Ensures the equality of all top level dataframes except dfs_*
    :return bool: True if equal
    """
    for var in self.__class__.instance_variables:
      expression1 = "self.%s is None" % (var)
      expression2 = "provider.%s is None" % (var)
      if eval(expression1) or eval(expression2):
        if eval(expression1) and eval(expression2):
          next
        else:
          return False
      else:
        expression = "isinstance(self.%s, pd.DataFrame)" % var
        if eval(expression):
          expression = "self.%s.equals(provider.%s)" % (var, var)
          if not eval(expression):
            return False
    return True

  def _makeGoTerms(self):
    df = self._makeDFFromCSV(FILENAME_GO_TERMS, is_index_geneid=True)
    df[cn.GO_TERM] = [v.strip() for v in df[cn.GO_TERM]]
    df = df.sort_values(cn.GO_TERM)
    return df

  @property
  def dfs_adjusted_read_count(self):
    """
    Creates the dataframe replicats adjusted for
    library size and gene length.
    :return list-pd.DataFrame:
        columns are times
    """
    if self._dfs_adjusted_read_count is None:
      self._dfs_adjusted_read_count =   \
          [self.normalizeReadsDF(df) for df in self.dfs_read_count]
    return self._dfs_adjusted_read_count

  @property
  def dfs_adjusted_read_count_wrt0(self):
    """
    Creates the dataframe replicats adjusted for
    library size and gene length normalized w.r.t. time 0.
    :return list-pd.DataFrame:
        columns are times
    """
    if self._dfs_adjusted_read_count_wrt0 is None:
      dfs = []
      for df in self.dfs_adjusted_read_count:
        sers = []
        for column in df.columns:
          sers.append(df[column] / df[T0])
        df_adj = pd.DataFrame(sers)
        df_adj.index = df.columns
        dfs.append(df_adj.T)
      self._dfs_adjusted_read_count_wrt0 = dfs
      #self._dfs_adjusted_read_count_wrt0 = [ df/df[T0] for df in self.dfs_adjusted_read_count]
    return self._dfs_adjusted_read_count_wrt0

  @property
  def dfs_centered_adjusted_read_count(self):
    def center(df):
      """
      Centers the dataframe along the rows
      :param pd.DataFrame:
      :return pd.DataFrame: same index, column
      """
      df_result = df.T
      df_result = df_result - df_result.mean()
      return df_result.T
    #
    if self._dfs_centered_adjusted_read_count is None:
      self._dfs_centered_adjusted_read_count =  \
          [center(df) for df in self.dfs_adjusted_read_count]
    return self._dfs_centered_adjusted_read_count

  def do(self, data_dir=cn.DATA_DIR):
    """
    Assigns values to the instance data.
    """
    persister = Persister(cn.DATA_PROVIDER_PERSISTER_PATH)
    if persister.isExist():
      provider = persister.get()
      self._setValues(provider=provider)
    else:
      # Gene categorizations
      self.df_ec_terms =  \
          self._makeDFFromCSV(FILENAME_EC_TERMS, is_index_geneid=True)
      self.df_ko_terms =  \
          self._makeDFFromCSV(FILENAME_KO_TERMS, is_index_geneid=True)
      self.df_kegg_pathways =  \
          self._makeDFFromCSV(FILENAME_KEGG_PATHWAYS,
          is_index_geneid=False)
      self.df_kegg_gene_pathways =  \
          self._makeDFFromCSV(FILENAME_KEGG_GENE_PATHWAY,
          is_index_geneid=True)
      # GO Terms
      self.df_go_terms = self._makeGoTerms()
      # Gene expression for state
      self.df_gene_expression_state = self._makeDFFromCSV(
          FILENAME_GENE_EXPRESSION_STATE, is_index_geneid=True)
      # Gene description
      self.df_gene_description = self._makeGeneDescriptionDF()
      # Stages matrix
      self.df_stage_matrix = self._makeStageMatrixDF()
      # Normalized data values
      self.df_normalized = self._makeNormalizedDF()
      # Raw readcounts
      self.dfs_read_count = self._makeReadCountDFS()
      # Hypoxia data
      self.df_hypoxia = self._makeHypoxiaDF()
      # Create mean and std dataframes
      self.df_mean = self._makeMeanDF()
      self.df_std = self._makeStdDF()
      self.df_cv = 100 * self.df_std / self.df_mean
      persister.set(self)
