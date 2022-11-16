#!/bin/bash
# Sets the python path to run codes in this repo
# Handle the sub-modules
PYTHONPATH=""
for d in common_python
  do
    cd $d
    PYTHONPATH=${PYTHONPATH}:`pwd`
    cd ..
  done
#
pushd xstate
cd python
PYTHONPATH=${PYTHONPATH}:`pwd`
export PYTHONPATH
popd
#
conda config --set changeps1 True 
source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate xstate
