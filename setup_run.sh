#!/bin/bash
# Sets the python path to run codes in this repo
# Handle the sub-modules
for d in common_python
  do
    cd $d
    source setup_run.sh
    cd ..
  done
#
pushd xstate
cd python
PYTHONPATH=${PYTHONPATH}:`pwd`
export PYTHONPATH
popd
