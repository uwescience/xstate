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
  dfs_adjusted_read_count_wrtT0 - adjusted w.r.t. libary size, gene length, time 0
     index: cn.GENE_ID
     column: time
  dfs_adjusted_read_count_wrtT0_log2 - in log2 units
     index: cn.GENE_ID
     column: time
  dfs_centered_adjusted_read_count - centers w.r.t. mean value of gene
     index: cn.GENE_ID
     column: time
  df_trn_unsigned - transcription regulation network
     column: cn.TF, cn.GENE_ID, cn.SIGN (always 1)
  df_trn_signed - TRN with sign for activation/inhibition
     column: cn.TF, cn.GENE_ID, cn.SIGN (1, -1)
  tfs - list of transcription factors

For unreplicated data, 
time columns have headings of the form T0, T1, ...
For replicated data, the format is "T%d.%d". The
latter is called the time replicated format.
"""

import common.constants as cn
import common_python.constants as cpn
from common_python.util.persister import Persister

import copy
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
FILENAME_TRN_UNSIGNED = "mtb.chip.network.drem"
FILENAME_TRN_SIGNED = "MTB-Signed-TRN"
NUM_REPL = 3
T0 = 0
MIN_LOG2_VALUE = -10
MILLION = 1e6
KILOBASE = 1e3  # Thousand bases
SEPARATOR = "."  # Separator used in time columns to indicate replication
# Columns
# Columns of files for transcription regulator networks
TRN_COLUMNS = [cn.TF, cn.GENE_ID, cn.SIGN]


##################### FUNCTIONS #####################
def calcRefDefault(df):
  return df[cn.TIME_0]


############################### CLASSES #########################
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
    "df_trn_signed",
    "df_trn_unsigned",
    "dfs_read_count",
    "_dfs_adjusted_read_count",
    "_dfs_adjusted_read_count_wrtT0",
    "_dfs_adjusted_read_count_wrtT0_log2",
    "_dfs_centered_adjusted_read_count",
    "tfs",
    ]

  def __init__(self, data_dir=cn.DATA_DIR,
      is_normalized_wrtT0=True,
      is_only_qgenes=True,
      is_display_errors=True,
      calcRef=calcRefDefault):
    """
    :param bool is_normalized_wrtT0: normalize data w.r.t. T0
        Otherwise, standardize values using the mean.
    :param bool is_only_qgenes: only include genes included in multi-hypothesis test
    :param calcRef Function: Calculate reference value for gene expression
        input: DataFrame
           columns: instances
           index: genes
           values: log2
        output: series
           index: genes
           values: log2
    """
    #
    if not is_normalized_wrtT0:
      raise ValueError("No longer support is_normalized_wrtT0 == False.")
    self.calcRef = calcRef
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
        if var in dir(provider):
          stmt = "self.%s = provider.%s" % (var, var)
        else:
          stmt = "self.%s = None" % var
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

  def removeReplicaStrings(self, dfs):
    """
    Removes the replica strings from columns of the dataframes.
    :param list-pd.DataFrame: columns are times
    :param list-pd.DataFrame:
    """
    new_dfs = copy.deepcopy(dfs)
    for df in new_dfs:
      new_columns = []
      for column in df.columns:
        splits = column.split(SEPARATOR)
        new_columns.append(splits[0])
      df.columns = new_columns
    return new_dfs

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
      dfs = self.removeReplicaStrings(self.dfs_centered_adjusted_read_count)
      dfs_new = [df.applymap(predicate) for df in dfs]
      return sum (dfs_new) / len(dfs)

  def _makeStdDF(self):
      """
      Creates a dataframe for the standard deviations.
      :return pd.DataFrame:
      """
      num_repl = self._getNumRepl()
      df_mean = self._makeMeanDF()
      dfs = self.removeReplicaStrings(
          self.dfs_centered_adjusted_read_count)
      df_std = (sum([dfs[n]*dfs[n]
          for n in range(num_repl)])  \
          - num_repl * df_mean * df_mean) /  \
          (num_repl - 1)
      return df_std.pow(1./2)

  def _reduceDF(self, df):
    """
    Reduces a dataframe indexed by GENE_ID to those 
    genes available in df_gene_description
    :param pd.DataFrame df: indexed by GENE_ID
    :return pd.DataFrame:
    """
    return pd.DataFrame([r for idx, r in df.iterrows()
        if idx in self.df_gene_description.index])

  def _makeReadCountDFS(self):
    """
    Creates a list of dataframes for each replication of the counts.
    :return list-pd.DataFrame:
      indexed by GENE_ID
      column names are in the time replicated format ("T%d.%d)
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
    int_map = {'A': 0, 'B': 1, 'C': 2}
    for repl in range(NUM_REPL):
      # Select the columns for this replication
      col_suffix = "_%s" % name_map[repl]
      column_names = [c for c in df.columns 
          if col_suffix in c]
      # Transform names into numbers
      new_names = {}
      for name in column_names:
        split_name = name.split("_")
        replica_letter = split_name[3]
        new_name = split_name[2][1:]
        suffix = int_map[replica_letter]
        new_names[name] = self._makeTime(
            int(new_name), suffix=suffix)
      df_repl = df[column_names]
      df_repl = df_repl.rename(columns=new_names)
      df_repl.index = df.index
      dfs.append(df_repl)
    #
    return dfs

  def calcRefPooled(self, df):
    """
    Calculates a reference value that is pooled value of
    normalized but not log2 counts.

    Parameters
    ----------
    df: DataFrame
        columns: times
        rows: genes
        values: log2 of counts
    
    Returns
    -------
    Series
        rows: genes
        values: log2 of counts
    """
    #df_denormalized = df.applymap(lambda v: 2**v)
    #ser = df_denormalized.mean(axis=1)
    #ser = ser.apply(lambda v: np.log2(v))
    ser = df.mean(axis=1)
    return ser

  def _getLog2NormalizedReadcounts(self):
    """
    Constructs a dataframe of log2 normalized reads.
    
    Returns
    -------
    DataFrame
        columns: timeponts
        index: genes
        values: log2 of counts
    """
    df = self._makeDFFromCSV(FILENAME_NORMALIZED)
    return df.set_index(cn.GENE_ID)

  def _makeNormalizedDF(self):
    """
    Transformation of the "normalized Read Counts" processed by DESeq2.
    Standardized the values for each gene.
    Drops rows where all columns are minimum values.
    Assumes that self.df_gene_expression_state has been initialized.
    Only includes genes that are expressed.
    :return pd.DataFrame:
        rows: gene
        columns: time
    """
    def defaultCalcRef(df):
        return df[cn.TIME_0]
    #
    df = self._getLog2NormalizedReadcounts()
    # Normalize w.r.t. the counts
    drops = []  # Rows to drop
    ser_ref = self.calcRef(df)
    for idx in df.index:
      values = df.loc[idx, :] - ser_ref.loc[idx]
      df.loc[idx, :] = [max(MIN_LOG2_VALUE,  v) for v in values]
      if all([v <= MIN_LOG2_VALUE for v in df.loc[idx, :]]):
        drops.append(idx)
    df = df.drop(index=drops) # Drop the 0 rows
    # Find genes to keep
    if self._is_only_qgenes:
      keep_genes = self.df_gene_expression_state.index
      df = df[df.index.isin(keep_genes)]
    #
    return df
 
  def normalizeReadsDF(self, df, suffix=None,
      is_time_columns=True):
    """
    Transforms read counts into features and units
    used in analysis.
    :param pd.DataFrame df:
    :param str suffix:  suffix appended to column names
      1. Adjusts read counts for each gene based on length and library size.
         Uses the metric RPKM - reads per kilobase millions
      2. Deletes unused features
    :param pd.DataFrame df: 
        index are genes, 
        columns are instances with the same values as self.df_normalized
        values are read counts
    :param bool is_time_columns: has columns of times
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
    if is_time_columns:
      df_result.columns = self.makeTimes(suffix=suffix)
    return df_result

  def _makeTime(self, int_time, suffix=None):
    result = "T%d" % int_time
    if suffix is not None:
      result = "%s%s%s" % (result, SEPARATOR, suffix)
    return result
  
  def makeTimes(self, suffix=None):
    """
    Creates the names for time.
    :return list-str:
    """
    if self.df_normalized is None:
      df = self._makeDFFromCSV(FILENAME_NORMALIZED)
    else:
      df = self.df_normalized
    return [self._makeTime(n, suffix=suffix) for n
        in range(len(df.columns))]

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
    df = self._makeDFFromCSV(FILENAME_GO_TERMS,
        is_index_geneid=True)
    df[cn.GO_TERM] = [v.strip() for v in df[cn.GO_TERM]]
    # Concatenate occurrences of GO terms
    dff = df.reset_index()
    dfg = dff.groupby(cn.GENE_ID)
    group_dct = dict(dfg.groups)
    term_dct = {cn.GENE_ID: [], cn.GO_TERM: []}
    for gene_id, indices in group_dct.items():
      terms = []
      for term in dff.loc[indices, cn.GO_TERM]:
        terms.append(term)
      term_dct[cn.GENE_ID].append(gene_id)
      term_dct[cn.GO_TERM].append("---".join(terms))
    df_result = pd.DataFrame(term_dct)
    return df_result

  @property
  def dfs_adjusted_read_count(self):
    """
    Creates the dataframe replicas adjusted for
    library size and gene length.
    :return list-pd.DataFrame:
        columns are times
    """
    if self._dfs_adjusted_read_count is None:
      self._dfs_adjusted_read_count = []
      for idx, df in enumerate(self.dfs_read_count):
        self._dfs_adjusted_read_count.append(
            self.normalizeReadsDF(df, suffix=idx))
    return self._dfs_adjusted_read_count

  @property
  def dfs_adjusted_read_count_wrtT0(self):
    """
    Creates the dataframe replicas adjusted for
    library size and gene length normalized w.r.t. time 0.
    :return list-pd.DataFrame:
        columns are times
    """
    if self._dfs_adjusted_read_count_wrtT0 is None:
      dfs = []
      for idx, df in  \
          enumerate(self.dfs_adjusted_read_count):
        sers = []
        t0_column = self.getT0s(df.columns)[0]
        for column in df.columns:
          sers.append(df[column] / df[t0_column])
        df_adj = pd.DataFrame(sers)
        df_adj.index = self.makeTimes(suffix=idx)
        dfs.append(df_adj.T)
      self._dfs_adjusted_read_count_wrtT0 = dfs
    return self._dfs_adjusted_read_count_wrtT0

  def getT0s(self, values):
    return [c for c in values if cn.TIME_0 in c]

  @property
  def dfs_adjusted_read_count_wrtT0_log2(self):
    """
    dfs_adjusted_read_count_wrtT0 in log2 units
    :return list-pd.DataFrame:
        columns are times, rows are genes
    """
    if self._dfs_adjusted_read_count_wrtT0_log2 is None:
      self._dfs_adjusted_read_count_wrtT0_log2 =  \
          [df.applymap(lambda v: np.log2(v) if v > 0 else 0)
          for df in self.dfs_adjusted_read_count_wrtT0]
    return self._dfs_adjusted_read_count_wrtT0_log2

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

  def getStageNames(self, ser_y):
    """
    Provides the list of stage names that correspond to the state indexes.

    Parameters
    ----------
    ser_y: Series
    
    Returns
    -------
    list-str
    """
    names = np.repeat(None, ser_y.max() + 1)
    last_name = None
    for timepoint, value in ser_y.iteritems():
      if "." in timepoint:
        new_timepoint = timepoint[0:-2]
      else:
        new_timepoint = timepoint
      names[value] = self.getStages(new_timepoint)
    return names

  def getStages(self, timepoints):
    """
    Returns the names of stages for the timepoints or integer time offsets.

    Parameters
    ----------
    timepoints: str/list-str/int/list-int
    
    Returns
    -------
    str/np.array-str
    """
    def convertIntToTimepoint(val):
      return "T%d" % val
    #
    if isinstance(timepoints, str):
      # str
      new_timepoints = [timepoints]
    elif isinstance(timepoints, int):
      # int
      new_timepoints = [convertIntToTimepoint(timepoints)]
    elif isinstance(timepoints[0], str):
      # list-str
      new_timepoints = timepoints
    elif isinstance(timepoints[0], int):
      new_timepoints = [convertIntToTimepoint(t) for t in timepoints]
    else:
      raise RuntimeError("Invalid type")
    #
    result = self.df_stage_matrix["name"].loc[new_timepoints].values
    return result

  def do(self, data_dir=cn.DATA_DIR):
    """
    Assigns values to the instance data.
    """
    persister = Persister(cn.DATA_PROVIDER_PERSISTER_PATH)
    done = False
    if persister.isExist():
      try:
        provider = persister.get()
        if "calcRef" in dir(provider):
          if str(self.calcRef) == str(provider.calcRef):
            self._setValues(provider=provider)
            done = True
      except AttributeError:
        pass
    if not done:
      # Gene categorizations
      self.df_ec_terms =  \
          self._makeDFFromCSV(FILENAME_EC_TERMS,
          is_index_geneid=True)
      self.df_ko_terms =  \
          self._makeDFFromCSV(FILENAME_KO_TERMS, 
          is_index_geneid=True)
      self.df_kegg_pathways =  \
          self._makeDFFromCSV(FILENAME_KEGG_PATHWAYS,
          is_index_geneid=False)
      self.df_kegg_gene_pathways =  \
          self._makeDFFromCSV(FILENAME_KEGG_GENE_PATHWAY,
          is_index_geneid=True)
      # Transcription Regulation Network
      self.df_trn_unsigned = self._makeDFFromCSV(
          FILENAME_TRN_UNSIGNED)
      self.df_trn_unsigned.columns = TRN_COLUMNS
      self.df_trn_signed = self._makeDFFromCSV(
          FILENAME_TRN_SIGNED)
      self.df_trn_signed.columns = TRN_COLUMNS
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
      # Transcription factors
      self.tfs = self.df_trn_unsigned[cn.TF].unique()
      self.tfs = list(set(self.tfs).intersection(
          self.dfs_adjusted_read_count[0].index))
      persister.set(self)
