import common.constants as xcn
from classifier import util_classifier
from common.trinary_data import TrinaryData
from common.data_provider import DataProvider
from common_python.classifier.feature_set import FeatureSet

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False


class TestFunctions(unittest.TestCase):

  def setUp(self):
    self.trinary = TrinaryData()
    self.provider = DataProvider()
    self.provider.do()

  def testCountTerms(self):
    if IGNORE_TEST:
      return
    TERM = "DNA replication"
    EXPECTED_COUNT = 2
    def test(terms, expected_count, fset=None):
      df_gene = self.provider.df_go_terms
      if fset is None:
        feature_set = FeatureSet(
            df_gene[xcn.GENE_ID][1:3])
        fset = FeatureSet(df_gene[xcn.GENE_ID][1:3])
      count = util_classifier.countTerms(fset, terms)
      self.assertEqual(count, expected_count)
    #
    test([TERM], EXPECTED_COUNT)
    test([TERM, TERM], 2*EXPECTED_COUNT)
    test(["DUMMY"], 0)
    #
    fset = FeatureSet(['Rv0981--Rv1332--Rv1828'])
    test(["DUMMY"], 0, fset=fset)

  def testCountTerms2(self):
    if IGNORE_TEST:
      return
    TERMS = ["a"]
    fset = FeatureSet(["Rv2009"])
    count1 = util_classifier.countTerms(fset, TERMS,
        is_include_module=False)
    count2 = util_classifier.countTerms(fset, TERMS,
        is_include_module=True)
    self.assertGreater(count2, count1)

  def testExtractAggregatedGene(self):
    if IGNORE_TEST:
      return
    GENES = ['Rv0981', 'Rv1332', 'Rv1828']
    AGGREGATED_GENE = xcn.GENE_SEPARATOR.join(GENES)
    genes = util_classifier.extractAggregatedGene(
        AGGREGATED_GENE)
    diff = set(GENES).symmetric_difference(genes)
    self.assertEqual(len(diff), 0)
    

if __name__ == '__main__':
  unittest.main()
