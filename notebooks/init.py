"""Initializes for using SBMLLint."""

import os
import sys


BASE_DIR = os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))
MINION_CODE_DIR = BASE_DIR
for directory in ["xstate", "python"]:
  MINION_CODE_DIR = os.path.join(MINION_CODE_DIR, directory)
sys.path.insert(0, MINION_CODE_DIR)
COMMON_CODE_DIR = os.path.join(BASE_DIR, "common_python")
sys.path.insert(0, COMMON_CODE_DIR)
