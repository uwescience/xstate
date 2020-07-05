"""Utilities for xstate classifiers."""

from common import constants as xcn
from common.data_provider import DataProvider
from common_python import constants as ccn


def extractAggregatedGene(gene):
  return gene.split(xcn.GENE_SEPARATOR)

def countTerms(feature_set, terms):
  """ 
  Counts the occurrences of terms in the GO terms
  of genes in the FeatureSet.

  Parameters
  ----------
  feature_set: FeatureSet
  terms: list-str

  Returns
  -------
  int
  """
  provider = DataProvider()
  provider.do()
  # Extract the genes
  genes = []
  [genes.extend(extractAggregatedGene(c)) for
      c in feature_set.list]
  # Compile the string of go terms for the genes
  df = provider.df_go_terms
  go_terms = [df[df[xcn.GENE_ID]==g][xcn.GO_TERM].values[
      0] for g in genes]
  go_str = "****".join(go_terms)
  count = sum([go_str.count(t) for t in terms])
  return count
