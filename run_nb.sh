#!/bin/csh
# Converts and executes an ipython notebook
# Notes: (1) File should not contain spaces
NOTEBOOK="$1.ipynb"
PYTHON="$1.py"
TMP1="/tmp/run_nb1.py"
if test -f "${NOTEBOOK}"; then
    echo "Creating ${PYTHON}$"
else
    echo "**Error. Cannot find ${NOTEBOOK}"
    exit -1
fi
#
jupyter nbconvert --to script "${NOTEBOOK}"
sed 's/^# In\[.*$/print("&")/' < ${PYTHON}  >  ${TMP1}
sed 's/^get_ipython().run/#&/' < ${TMP1}  > ${PYTHON}
python ${PYTHON}
