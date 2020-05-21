#!/bin/bash
# Should be run from this directory
# Runs processing in parallel for each state
cd ../../../
echo `pwd`
source setup_run.sh
cd xstate/python
for s in 0 1 2 3 4 5
  do
    echo "***Running class $s"
    python classifier/main_feature_equivalence_calculator.py $s --restart &
  done
