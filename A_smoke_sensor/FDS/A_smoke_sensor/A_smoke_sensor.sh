#!/bin/sh

# Execution of fire simulation and python script to supply the Smoke Sensor


#===================CHANGE PARAMETERS AS NECEASSARY=============================
chids=(
'smoke_sensor'
)

quantities=(
'EXTINCTION COEFFICIENT'
)

locations=(
'Z',2.25
)
#==================AUTOMATIC PART - NO CHANGES NECEASSARY=======================

# cd ..
# for chid in "${chids[@]}"
# do
#   echo Running fds simulation...
#   fds $chid.fds
#   wait
# done

cd A_smoke_sensor

for location in "${locations[@]}"
do
  IFS=","; set $location;
  # $1 = slice_dim as str and $2 = slice_coord as float

  for quantity in "${quantities[@]}"
  do
     echo Slicefile quantity: "$quantity"
     echo Slicefile dimension: "$1"
     echo Slicefile coordinate: $2

     python3 0_main_script.py -q "$quantity" -d "$1" -c $2

  done
done
