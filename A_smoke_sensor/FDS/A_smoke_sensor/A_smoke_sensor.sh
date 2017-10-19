#!/bin/sh

# Execution of the python scripts to supply the Smoke Sensor


#===================CHANGE PARAMETERS AS NECEASSARY=============================
quantities=(
'SOOT OPTICAL DENSITY'
)

locations=(
'Z',2.25
)
#===============================================================================


#### AUTOMATIC PART - NO CHANGES NECEASSARY ####

for location in "${locations[@]}"
do
  IFS=","; set $location;
  # $1 = slice_dim as str and $2 = slice_coord as float

  for quantity in "${quantities[@]}"
  do
     echo Slicefile quantity: "$quantity"
     echo Slicefile dimension: "$1"
     echo Slicefile coordinate: $2

     python 0_slice2ascii.py -q "$quantity" -d "$1" -c $2
     wait
     python 1_meshgrid.py
     wait
     python 2_consolidate_meshes.py
     wait
     python 3_smoke_factor_grids.py

  done
done
