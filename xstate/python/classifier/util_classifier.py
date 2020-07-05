"""Utilities for xstate classifiers."""

from common import constants as xcn
from common.data_provider import DataProvider
from common_python import constants as ccn


def extractAggregatedGene(gene):
  return gene.split(xcn.GENE_SEPARATOR)

def countTerms(fset, terms,
     is_include_module=True):
  """ 
  Counts the occurrences of terms in the GO terms
  of genes in the FeatureSet.

  Parameters
  ----------
  fset: FeatureSet
  terms: list-str
  is_include_module: bool
      consider all genes in modules activated by a
      gene in fset

  Returns
  -------
  int
  """
  provider = DataProvider()
  provider.do()
  # Extract the genes
  genes = []
  [genes.extend(extractAggregatedGene(c)) for
      c in fset.list]
  if is_include_module:
    new_genes = []
    tfs = list(set(
        provider.df_trn_unsigned[xcn.TF].tolist()))
    for gene in genes:
      if gene in tfs:
        sel = provider.df_trn_unsigned[xcn.TF] == gene
        df = provider.df_trn_unsigned[sel]
        new_genes.extend(df[xcn.GENE_ID].tolist())
    genes.extend(new_genes)
    genes = list(set(genes))
  # Compile the string of go terms for the genes
  df = provider.df_go_terms
  indices = [df[df[xcn.GENE_ID]==g].index.tolist()
      for g in genes]
  indices = [t for l in indices for t in l]
  go_terms = df.loc[indices, xcn.GO_TERM].to_list()
  go_str = "****".join(go_terms)
  count = sum([go_str.count(t) for t in terms])
  return count
