#!/bin/sh

# Execution of the python scripts to supply the Smoke Sensor

python 0_slice2ascii.py
wait
python 1_meshgrid.py
wait
python 2_consolidate_meshes.py
wait
python 3_process_collect.py
wait
python 4_sfgrids.py
