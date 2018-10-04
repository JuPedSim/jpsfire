#!/bin/sh

# Execution of fire simulation and python script to supply the Smoke Sensor


#===================CHANGE PARAMETERS AS NECEASSARY=============================
chids=(
'smoke_sensor'
)

quantities=(
'CO2'
'CO'
'HCN'
'TEMPERATURE'
'EXTINCTION COEFFICIENT'
)

locations=(
'Z',2.25
)
#==================AUTOMATIC PART - NO CHANGES NECEASSARY=======================

for quantity in "${quantities[@]}"
do
  echo Slicefile quantity: "$quantity"
  python3 preprocess_fds.py -q "$quantity" -p
done
