#!/bin/bash
# Create cross validation data. Run from this directory.
cd ..
python tools/cross_validation_data.py --clean=True
python tools/cross_validation_data.py --start=0 --end=2 &
python tools/cross_validation_data.py --start=2 --end=4 &
python tools/cross_validation_data.py --start=4 --end=6 &
python tools/cross_validation_data.py --start=6 --end=8 &
