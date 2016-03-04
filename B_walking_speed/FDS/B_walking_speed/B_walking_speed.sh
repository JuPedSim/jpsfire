#!/bin/sh

# Execution of the python scripts to supply the Smoke Sensor

quantities=(
'EXTINCTION COEFFICIENT'
)

# Please specify the coordinate as float!
locations=(
'Z',1.8
'Z',4.8
)

for location in "${locations[@]}"
do
  IFS=","; set $location;
  # $1 = slice_dim as str and $2 = slice_coord as float

  for quantity in "${quantities[@]}"
  do
     echo Slicefile quantity: "$quantity"
     echo Slicefile dimension: "$1"
     echo Slicefile coordinate: $2

     python 0_slice2ascii.py "$quantity" "$1" $2
     wait
     python 1_meshgrid.py
     wait
     python 2_consolidate_extinction_grids.py

  done
done
