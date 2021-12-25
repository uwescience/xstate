import sys
import os
this_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_dir)
#
import add_path
import common.constants as cn
from common_python.util import dataframe
from common_python.util.persister import Persister
from tools import make_classification_data

import os
import shutil
import unittest

IGNORE_TEST = False
IS_PLOT = False
IS_CHANGED = False # ClassificationData namespace has changed
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE_PATH = os.path.join(TEST_DIR, make_classification_data.DATA_FILE)
DATA_FILE_PATH_TEST = os.path.join(
    TEST_DIR, "make_classification_data_save.pcl")
PERSISTER = Persister(DATA_FILE_PATH_TEST)
FILES = [DATA_FILE_PATH]
data= make_classification_data.ClassificationData(
    persister_path=DATA_FILE_PATH_TEST)
# Creates a file with the desired initializtions, if necessary
# Note, if there are changes in the state variables, you
# must set IS_CHANGED to True.
if not PERSISTER.isExist() or IS_CHANGED:
  data.initialize()
  data.serialize()


class TestClassificationData(unittest.TestCase):

  def deleteFiles(self):
    for ffile in FILES:
      if os.path.isfile(ffile):
        os.remove(ffile)

  def setUp(self):
    self.deleteFiles()
    shutil.copyfile(DATA_FILE_PATH_TEST, DATA_FILE_PATH)
    self.data = make_classification_data.ClassificationData(
        persister_path=DATA_FILE_PATH)

  def _isSameDict(self, dct1, dct2):
    self.assertEqual(len(dct1.keys()), len(dct2.keys()))
    for key, value in dct1.items():
      if "equals" in dir(value):
        self.assertEqual(value.equals(dct2[key]))
      else:
        self.assertEqual(value, dct2[key])

  def tearDown(self):
    self.deleteFiles()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertIsNotNone(self.data.persister)

  def testSetNamespace(self):
    if IGNORE_TEST:
      return
    NAMESPACE = {"a": 1, "b": "hello"}
    for key, value in NAMESPACE.items():
      self.data.namespace_dct[key] = value
    newNamespace = {}
    self.data.setNamespace(newNamespace)
    self._isSameDict(NAMESPACE, newNamespace)
 
  def testInitialize(self):
    if IGNORE_TEST:
      return
    self.data.persister.remove()
    self.data.initialize()
    self.data.setNamespace(globals())
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
 
  def testCalcAccuracy(self):
    if IGNORE_TEST:
      return
    self.data.deserialize()
    df = self.data.calcAccuracy(is_debug=True)
    diff = set(df.columns).symmetric_difference(
        make_classification_data.DF_ACCURACY_COLUMNS)
    self.assertEqual(len(diff), 0)
    self.assertGreater(len(df), 0)
    
      

if __name__ == '__main__':
  unittest.main()
