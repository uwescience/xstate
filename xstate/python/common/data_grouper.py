'''Groups Data According to Trinary Values'''

from common import constants as cn
from common.data_provider import DataProvider
from common import transform_data

import pandas as pd
import numpy as np


class DataGrouper(object):

  def __init__(self, df_trinary=None):
    """
    :param pd.DataFrame df: trinary valued DF (has values -1, 0, 1)
    """
    if df_trinary is None:
      provider = DataProvider()
      provider.do()
      df_trinary = transform_data.makeTrinaryData(
          is_include_nan=False)
    self.df_trinary = df_trinary
    self.df_group = None  # Dataframe describing groups
    self.df_gene_group = None  # Genes by group

  def _makeGroupDF(self):
    """
    Forms groups of rows with identical discrete values.
    :param pd.DataFrame df:
    :return updates self.df_group:
       index - trinary differential expressions
           tuple of {-1,0, 1}
       GENE_IDS - list of GENE_ID
       CNT_GROUP - count of members in group
       CNT_UP - count of columns with 1 in group
       CNT_DOWN - count of columns with -1 in group
    """
    dfg = self.df_trinary.groupby(self.df_trinary.columns.tolist())
    groups = {k: list(v) for k, v in dfg.groups.items()}
    expressions = list(groups.keys())
    gene_ids = [v.tolist() for v in dfg.groups.values()]
    cnt_groups = [len(v) for v in groups.values()]
    cnt_ups = [k.count(1) for k in groups.keys()]
    cnt_downs = [k.count(-1) for k in groups.keys()]
    self.df_group = pd.DataFrame({
        cn.CNT_GROUP: cnt_groups,
        cn.CNT_UP: cnt_ups,
        cn.CNT_DOWN: cnt_downs,
        cn.GENE_IDS: gene_ids,
        }, index=expressions)

  def _makeGeneGroupDF(self, is_include_0_group=False):
    """
    Creates a dataframe with gene groups.
    :param bool is_include_0_group: include genes without DE
    Constructs self.df_gene_group.
        index: GENE_ID
        Columns: GROUP
    """
    if self.df_group is None:
      self._makeGroupDF()
    dfs = []
    if is_include_0_group:
      exclude_group = ()
    else:
      a_group = self.df_group.index[0]
      exclude_group = tuple(np.repeat(0, len(a_group)))
    for group, row in self.df_group.iterrows():
      if group == exclude_group:
        continue
      gene_ids = row[cn.GENE_IDS]
      groups = [group] * len(gene_ids)
      df = pd.DataFrame({
        cn.GENE_ID: gene_ids,
        cn.GROUP: groups
        })
      dfs.append(df)
    self.df_gene_group = pd.concat(dfs)
    self.df_gene_group = self.df_gene_group.set_index(cn.GENE_ID)
      
  def _pruneSize(self, min_size=2, is_include_zeros=False):
    """
    Eliminates small groups and 0 occurrences.
    :param int min_size:
    """
    if not is_include_zeros:
      for idx, values in enumerate(self.df_group.index):
        if all([v == 0 for v in values]):
          self.df_group = self.df_group.drop(values, axis=0)
          break
    self.df_group = self.df_group[
        self.df_group[cn.CNT_GROUP] >= min_size]

  def do(self, **kwargs):
    """
    Prepares the grouper data frame
    """
    self._makeGroupDF()
    self._pruneSize(**kwargs)
    self._makeGeneGroupDF()
