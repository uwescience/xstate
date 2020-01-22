'''Creates and manipulates term matrices, where terms are columns'''

from common import constants as cn
from common.data_provider import DataProvider
from common.data_grouper import DataGrouper
from common_python.statistics import util_statistics
from plots import util_plots
from common_python.text import util_text

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
from scipy.spatial import distance
from scipy.cluster.hierarchy import linkage, fcluster

HIGH_SL = 1e6

class TermMatrix(object):

  """
  The core dataframe is the term matrix, self.df_matrix. Its columns
  are terms; the rows are groups of correlated genes. A group is
  a tuple of trinary values indicating when that terms is expressed.
  """

  def __init__(self, term_column=cn.GO_TERM, 
      is_plot=True, **kwargs):
    """
    :param str term_column: column in go_terms to use for text
    :param dict **kwargs: arguments to DataGrouper
    """
    self._is_plot = is_plot
    self.provider = DataProvider()
    self.provider.do()
    self.grouper = DataGrouper(**kwargs)
    self.grouper.do(min_size=1)
    self.df_matrix = self._makeMatrix(term_column)
    self.df_gene_term = self._makeGeneTerm()

  def _makeTermGroup(self, term_column=cn.GO_TERM):
    """
    :param str term_column: column in go_terms to use for text
    :return pd.DataFrame:
      index - group (time intervals with trinary values)
      column - Term
    """
    df = self.grouper.df_gene_group.merge(
        self.provider.df_go_terms, left_index=True, 
        right_index=True, how='inner')
    if term_column == cn.INDEX:
      df_term = df[[cn.GROUP]].copy()
      df_term[term_column] = df.index
    else:
      df_term = df[[cn.GROUP, term_column]].copy()
    df_term = df_term.set_index(cn.GROUP)
    return df_term

  def _makeMatrix(self, term_column=cn.GO_TERM):
    """
    :param str term_column: column in go_terms to use for text
    :return pd.DataFrame: matrix with the terms
    """
    df_term = self._makeTermGroup(term_column=term_column)
    df_result = util_text.makeTermMatrix(df_term[term_column])
    return df_result

  def _makeGeneTerm(self):
    """
    Finds the genes and terms that co-occur at the same times.
    :return pd.DataFrame:
        cn.GROUP - trinary values for genes at times
        cn.TERM - list of GO terms
        cn.GENE_D - list of genes
        cn.CNT_TERM - count of GO terms
        cn.CNT_GENE - count of genes
        cn.CNT_REGULATED - count of times up- down-regulated
    """
    def makeGroupedDF(df):
      df = df.reset_index()
      return df.groupby(cn.GROUP)
    def extract(df, key, col):
      return df.loc[[key], col].values.tolist()
    #
    df_term = self._makeTermGroup()
    df_gene = self.grouper.df_gene_group
    df_gene = df_gene.reset_index()
    df_gene = df_gene.set_index(cn.GROUP)
    dfg_term = makeGroupedDF(df_term)
    dfg_gene = makeGroupedDF(self.grouper.df_gene_group)
    # Find the keys in common
    keys_term = [k for k in dfg_term.groups]
    keys_gene = [k for k in dfg_gene.groups]
    keys_common = set(keys_term).intersection(keys_gene)
    dict_df = {cn.GROUP: [], cn.TERM: [], cn.GENE_ID: []}
    for key in keys_common:
      dict_df[cn.GROUP].append(key)
      dict_df[cn.TERM].append(extract(df_term, key, cn.GO_TERM))
      dict_df[cn.GENE_ID].append(extract(df_gene, key, cn.GENE_ID))
    df_result = pd.DataFrame(dict_df)
    df_result[cn.CNT_GENE] = [len(g) for g in df_result[cn.GENE_ID]]
    df_result[cn.CNT_TERM] =  [len(t) for t in df_result[cn.TERM]]
    df_result[cn.CNT_REGULATED] =   \
        [sum([np.abs(x) for x in g]) for g in df_result[cn.GROUP]]
    return df_result

  def makeAggregationMatrix(self, predicates):
    """
    Creates a matrix with columns the same as df_matrix
    and row i that is the summation of the values in rows
    that satisfy predicate i.
    :param list-BooleanFunc predicates: predicate on group tuples
    """
    columns = self.df_matrix.columns
    column_values = {c.strip(): [] for c in columns}
    for pos, predicate in enumerate(predicates):
      row = np.repeat(0.0, len(columns))
      row = row.reshape(1, len(columns))
      # TODO: Fix predicates
      if False:
        for group in self.df_matrix.index:
          if predicate(group):
            values = np.array(self.df_matrix.loc[[group], :])
            row += values
      for group in self.df_matrix.index:
        if group[pos] == 1:
          values = np.array(self.df_matrix.loc[[group], :])
          row += values
      row = row.reshape(len(columns))
      # Add the row for this predicate
      for idx, col in enumerate(columns):
        column_values[col].append(row[idx])
    return pd.DataFrame(column_values)

  # TODO: Fix use of predicates
  def plotAggregation(self, predicates, min_val=0, 
      is_include_ylabels=True):
    df = self.makeAggregationMatrix(predicates)
    df = df.applymap(lambda v: 0 if v < min_val else v)
    df = df.applymap(lambda v: np.nan if np.isclose(v, 0) else v)
    # Drop columns that are all nans
    for col in df.columns:
      if all([np.isnan(v) for v in df[col]]):
        del df[col]
    # Construct the plot
    plt.subplot(1, 2, 2)
    heatmap = plt.pcolor(df.transpose(), cmap='jet')
    if is_include_ylabels:
      ax = plt.gca()
      ax.set_yticks(np.arange(len(df.columns))+0.5)
      ax.set_yticklabels(df.columns, fontsize=8)
    plt.title("Term Counts")
    plt.colorbar(heatmap)
    plt.show()

  def makeTimeAggregationMatrix(self, is_up_regulated=True):
    """
    Creates a matrix with columns the same as df_matrix
    and row i that is the summation of the values in rows
    that satisfy predicate i.
    :param bool is_up_regulated:
    """
    if is_up_regulated:
      direction = 1
    else:
      direction = -1
    columns = self.df_matrix.columns
    column_values = {c.strip(): [] for c in columns}
    for time in range(cn.NUM_TIMES):
      row = np.repeat(0.0, len(columns))
      row = row.reshape(1, len(columns))
      for group in self.df_matrix.index:
        if group[time] == direction:
          values = np.array(self.df_matrix.loc[[group], :])
          row += values
      row = row.reshape(len(columns))
      # Add the row for this predicate
      for idx, col in enumerate(columns):
        column_values[col].append(row[idx])
    return pd.DataFrame(column_values)

  def calcClusters(self, max_distance=1, is_up_regulated=True):
    """
    Calculates log significance levels and clusters.
    :param float max_distance: maximum distance between clusters,
        otherwise merged
    :return pd.DataFrame, ndarray, pd.Series:
       df_log - log of significance level
       row_linkage - linkage matrix
           See https://stackoverflow.com/questions/9838861/scipy-linkage-format
       ser_cluster - cn.GROUP (indexed by term)
    """
    df = self.makeTimeAggregationMatrix(
        is_up_regulated=is_up_regulated)
    # Remove rows with zero variance
    df_filtered = util_statistics.filterZeroVarianceRows(df.T)
    # Compute significance levels
    df_log = util_statistics.calcLogSL(df_filtered, round_decimal=3)
    df_log = df_log.applymap(lambda v:
        HIGH_SL if np.isnan(v) else v)
    # Compute the clusters
    log_arrays = np.asarray(df_log)
    row_linkage = linkage(
        distance.pdist(log_arrays), method='average')
    ser_cluster = pd.Series(fcluster(row_linkage, 0.1, criterion="distance"))
    ser_cluster.index = df_log.index
    #
    return df_log, row_linkage, ser_cluster

  # Include state transitions
  # Note how clusters relate to state observations
  def plotTimeAggregation(self, is_up_regulated=True):
    """
    Plots aggregation of groups over time.
    :param bool is_include_ylabels:
    :param bool is_up_regulated:
    """
    df_log, row_linkage, ser_cluster =  \
        self.calcClusters(is_up_regulated=is_up_regulated)
    # Heatmap
    cg = sns.clustermap(df_log, row_linkage=row_linkage, 
        col_cluster=False, cbar_kws={"ticks":[0,5]}, cmap="Blues")
    # Construct a cluster map
    #cg = sns.clustermap(df_log, col_cluster=False, 
    #    cbar_kws={"ticks":[0,5]}, cmap="Blues")
    # Set the labels
    cg.ax_heatmap.set_xlabel("Time")
    if is_up_regulated:
      direction = "Up"
    else:
      direction = "Down"
    title = "-log10 zscores of %s-regulated term counts" % (direction)
    cg.ax_heatmap.set_title(title)
    xticks = cg.ax_heatmap.get_xticks() - 0.5  # Correct tick position
    cg.ax_heatmap.set_xticks(xticks)
    cg.ax_heatmap.set_yticks([])
    cg.ax_heatmap.set_yticklabels([])
    # Add the state transitions
    util_plots.plotStateTransitions(ymax=len(df_log),
        ax=cg.ax_heatmap, is_plot=False)
    #
    if self._is_plot:
      plt.show()
