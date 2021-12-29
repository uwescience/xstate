"""Creates Trinary Data for all Samples"""


import common.constants as cn
from common import trinary_data

import argparse
import os


def run():
  data = trinary_data.getSampleData()
  data.serialize()


if __name__ == '__main__':
  msg = """
Create trinary valued files for samples.
Samples are plaed in %s.
""" % cn.TRINARY_SAMPLES_DIR
  parser = argparse.ArgumentParser(description=msg)
  run()
