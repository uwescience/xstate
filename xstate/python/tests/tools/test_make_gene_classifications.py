"""Tests for make_gene_classifications."""

import common.constants as cn
from tools import make_gene_classifications
from common_python.testing import helpers

import os
import unittest

IGNORE_TEST = False


class TestMakeGeneClassifications(unittest.TestCase):

  def deleteFiles(self):
    for path in self.getPaths():
      if os.path.isfile(path):
        os.remove(path)

  def setUp(self):
    self.deleteFiles()

  def tearDown(self):
    self.deleteFiles()

  def getPaths(self):
    return [os.path.join(cn.TEST_DIR, f)
        for f in make_gene_classifications.FILES_KEGG]

  def testMake(self):
    if IGNORE_TEST:
      return
    make_gene_classifications.make(max_count=5, directory=cn.TEST_DIR)
    for path in self.getPaths():
      self.assertTrue(os.path.isfile(path))
      

if __name__ == '__main__':
  unittest.main()
