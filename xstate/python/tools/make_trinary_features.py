"""Creates Trinary Data for all Samples"""


import common.constants as cn
from common import trinary_data

import argparse
import os


def run():
  sample_data = trinary_data.getSampleData()
  trinary_data.mkFeatureMatrices(sample_data)


if __name__ == '__main__':
  msg = """
Create trinary valued files for samples.
Samples are plaed in %s.
""" % cn.TRINARY_SAMPLES_DIR
  parser = argparse.ArgumentParser(description=msg)
  run()
