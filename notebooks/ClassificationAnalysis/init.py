import os
import sys

base_dir = "/home/xstate"
xstate_dir = str(base_dir)
for directory in ["xstate", "python"]:
  xstate_dir = os.path.join(xstate_dir, directory)
sys.path.insert(0, xstate_dir)
common_dir = os.path.join(base_dir, "common_python")
sys.path.insert(0, common_dir)
