#!/bin/sh

# Execution of the python scripts to supply the Smoke Sensor

python 0_slice2ascii.py
wait
python 1_meshgrid.py
wait
python 2_consolidate_toxgrids.py
