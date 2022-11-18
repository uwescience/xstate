#!/bin/bash
source bin/alias.sh
DIR=`pwd`/xst
PATHS=$DIR:`pwd`/xstate/python:${DIR}/site-packages
source setup_run.sh ${PATHS}
export PYTHONPATH
source ${DIR}/bin/activate
