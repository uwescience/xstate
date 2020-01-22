"""Analyzes Terms in terms of the underlying gene structure and comparisons with other terms."""

"""
A term ontology is a classification of genes. Examples include: GO (gene ontology),
KO (KEGG Orthology), KEGG Pathway, and EC (Enzyme Commission). A term ontology
is a many-to-many relationship between genes and terms. A gene need not have
a corresponding term in a term ontology.
"""

from common import constants as cn
from common_python import constants as cpn
from common_python.plots import util_plots
from common.data_provider import DataProvider

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
from scipy.spatial import distance
from scipy.cluster.hierarchy import linkage, fcluster


class TermAnalyzer(object):

  def __init__(self, df_term, provider=None, is_plot=True):
    """
    :param pd.DataFrame df_term: a one column DataFrame and indexed by cn.GENE_ID
    :param bool is_plot:
    """
    self._is_plot = is_plot
    if provider is None:
      self.provider = DataProvider()
      self.provider.do()
    else:
      self.provider = provider
    self.df_term = df_term
    self.ontology = self._getOntology()

  def _getOntology(self):
    columns = [c for c in self.df_term.columns]
    return columns[0]

  def makeAnalyzerMatrix(self):
    """
    An analyzer matrix is a dataframe with columns that are terms (plus "Missing"),
    indexed by GENE_ID, and values are either a count or np.nan
    :return pd.DataFrame: analyzer matrix
    """
    # Create a matrix of expressed genes
    df_expressed = pd.DataFrame({
        cn.GENE_ID: self.provider.df_normalized.index,
        })
    df_expressed[self.ontology] = np.nan
    df_expressed[cn.COUNT] = np.nan
    # Matrix of terms
    df_term = self.df_term[self.df_term.index.isin(
        df_expressed[cn.GENE_ID])].copy()
    df_term = df_term.reset_index()
    df_term[cn.COUNT] = 1
    df_term = df_term.drop_duplicates()
    # Ensure all expressed genes are present
    gene_expressed = set(df_expressed[cn.GENE_ID].tolist())
    gene_term = set(df_term[cn.GENE_ID].tolist())
    gene_excluded = gene_expressed.difference(gene_term)
    df_expressed_excluded = df_expressed[
        df_expressed[cn.GENE_ID].isin(gene_excluded)].copy()
    df1_term = pd.concat([df_term, df_expressed_excluded])
    df_matrix = df1_term.pivot(index=cn.GENE_ID,
        columns=self.ontology, values=cn.COUNT)
    # Update the name of the nan column, if any
    columns = df_matrix.columns.tolist()
    try:
      idx = columns.index(np.nan)
    except ValueError:
      idx = -1
    if idx >= 0:
      columns[idx] = "Missing"
      df_matrix.columns = columns
    return df_matrix

  def plotTermHeatmap(self, **plot_opts):
    """
    Plots of heatmap of the expressed genes (x-axis) vs. terms (y-axis). Values
    are binary.
    """
    # Create a matrix of expressed genes
    df_matrix = self.makeAnalyzerMatrix()
    #df_matrix = df_matrix.fillna(0)
    opts = dict(plot_opts)
    opts = {cpn.PLT_CMAP: "Greys"}
    util_plots.plotCategoricalHeatmap(df_matrix, **opts)
