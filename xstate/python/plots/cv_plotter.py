"""Plots reated to Coefficient of Variation Analysis."""

from common import constants as cn
from common.data_provider import DataProvider

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


class CVPlotter():

  def __init__(self, provider=None, is_plot=True):
    if provider is None:
      self._provider = DataProvider()
      self._provider.do()
    else:
      self._provider = provider
    self._is_plot = is_plot

  def heatMap(self, min_cv=0):
    """
    Plots a heatmap of the coefficient of variations.
    :param pd.DataFrame df_cv: CVs
    :param float min_cv: minimum CV to consider
    """
    plt.figure(figsize=(16, 10))
    df = self._provider.df_cv
    # Rename columns to their hours
    ax = plt.gca()
    ax.set_xticks(np.arange(len(df.columns))+0.5)
    ax.set_xticklabels(df.columns)
    df = df.applymap(lambda v: v if v >= min_cv else np.nan)
    heatmap = plt.pcolor(df, cmap='jet')
    plt.colorbar(heatmap)
    plt.xlabel("times")
    plt.ylabel("gene")
    plt.title("Coefficient of Variation > %d percent" % min_cv)
    if self._is_plot:
      plt.show()

  def readsAndDO(self):
    """
    Plots the following lines for the hours of the experiments:
      Average CV of genes
      CV of dissolved oxygen (DO)
      Avg dissolved oxygen
    """
    hours = self._provider.df_hypoxia[cn.HOURS]
    means = self._provider.df_hypoxia[cn.MEAN]
    error_bars = [2*s for s in self._provider.df_hypoxia[cn.STD]]
    plt.errorbar(hours, means, yerr=error_bars, marker="o")
    ax = plt.gca()
    # Plot CVs of genes
    ser = self._provider.df_cv.mean()  # Average over geans
    ax.plot(hours, ser.values, linestyle='dashed', marker="o",
        color='r')
    plt.xlabel("hours")
    plt.ylabel("DO/CV")
    plt.legend(["CV for read counts", "DO +/- 2 std"])
    if self._is_plot:
      plt.show()
