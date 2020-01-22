from plots import cv_plotter
import common.constants as cn

import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False


class TestCVPlotter(unittest.TestCase):

  def setUp(self):
    self.plotter = cv_plotter.CVPlotter(is_plot=False)

  def testHeatMap(self):
    if IGNORE_TEST:
      return
    # Only smoke tests
    self.plotter.heatMap()

  def testReadsAndDO(self):
    # Only smoke tests
    self.plotter.readsAndDO()


if __name__ == '__main__':
  unittest.main()
