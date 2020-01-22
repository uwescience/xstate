<img src="https://api.travis-ci.org/uwescience/xstate.svg?branch=master" width="100"/>
   
# xstate: State Analysis for Gene Expression

## Developer Usage Notes
1. Execute ``source setup_env.sh`` at the directory root
1. Run python codes from the directory src/python

## Installation
The installation assumes that you have installed ``git``
and ``anaconda`` (or ``miniconda``).

- ``git clone --recurse-submodules https://github.com/uwescience/xstate.git``
- ``cd xstate``
- ``bash setup.sh``

Windows users should manually perform the steps in setup.sh.

## Usage

For all use cases:
1. Start a new terminal session
1. ``cd xstate``
1. ``conda activate xstate``

When done, ``conda deactivate xstate``


Update the classifiers:
1. From the top level repository directory, ``python xstate/xstate/python/tools/make_svm_classifier.py``. (On windows, use ``\`` instead of ``/``.

Run a classification:
1. Start a new terminal session, change directory to ``xstate``, and activate the virtual environment.
1. From the top level repository directory, ``python xstate/xstate/python/tools/classify_expression.py <data_file>``. (On windows, use ``\`` instead of ``/``. Note: It is assumed that \<data file\> is
in the ``samples`` directory under ``data``. To change the directory,
use the ``--dir`` option.
