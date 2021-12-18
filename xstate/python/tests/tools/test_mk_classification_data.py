import sys
import os
this_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_dir)
#
import add_path
import common.constants as cn
from tools import mk_classification_data

import os
import unittest

IGNORE_TEST = True
IS_PLOT = True
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE_PATH = os.path.join(TEST_DIR, mk_classification_data.DATA_FILE)
FILES = [DATA_FILE_PATH]


class TestClassificationData(unittest.TestCase):

  def deleteFiles(self):
    for ffile in FILES:
      if os.path.isfile(ffile):
        os.remove(ffile)

  def setUp(self):
    self.deleteFiles()
    self.data = mk_classification_data.ClassificationData(
        dir_path=DATA_FILE_PATH)

  def tearDown(self):
    self.deleteFiles()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertIsNotNone(self.data.persister)

  def testInitialize(self):
    # TESTING
    self.data.initialize()
    for name, value in self.data.namespace.items():
      globals()[name] = value
    import pdb; pdb.set_trace()
    for classifier_key, classifier_clf in CLASSIFIER_DCT.items():
        # Got expected genes in the feature vector?
        for gene_key, gene_list in GENE_DCT.items():
            if gene_key == classifier_key[1]:
                assert(len(gene_list) >= len(classifier_clf._df_X.columns))
        # Got expected data
        for data_key in DATA_DCT.keys():
            if data_key == classifier_key[0]:
                # Find common genes
                columns = GENE_DCT[classifier_key[1]]
                test_df = dataframe.subset(DATA_DCT[data_key].df_X, columns)
                newColumns = list(test_df.columns)
                newColumns.sort()
                test_df = test_df[newColumns]
                df_X = classifier_clf._df_X.copy()
                df_X = df_X[newColumns]
                self.assertTrue(test_df.equals(df_X))
      

if __name__ == '__main__':
  unittest.main()
