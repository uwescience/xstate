import sys
import os
this_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_dir)
#
import add_path
import common.constants as cn
from tools import classify_expression

import os
import unittest

IGNORE_TEST = False
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.join(TEST_DIR, "test_classify_expression.csv")
TEST_EXPRESSION_FILE = "AM_MDM_Mtb_transcripts_DEseq.csv"


class TestClassifyExpression(unittest.TestCase):

  def deleteFiles(self):
    if os.path.isfile(TEST_PATH):
      os.remove(TEST_PATH)

  def setUp(self):
    self.deleteFiles()

  def tearDown(self):
    self.deleteFiles()

  def testMake(self):
    if IGNORE_TEST:
      return
    classify_expression.processExpressionData(
        TEST_EXPRESSION_FILE, output_path=TEST_PATH,
        is_display_errors=False)
    self.assertTrue(os.path.isfile(TEST_PATH))
      

if __name__ == '__main__':
  unittest.main()
