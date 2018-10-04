#!/bin/sh

# Execution of fire simulation and python script to supply the Smoke Sensor


#===================CHANGE PARAMETERS AS NECEASSARY=============================
chids=(
'smoke_sensor'
)

quantities=('CO2','EXTINCTION COEFFICIENT', 'CO', 'HCN', 'TEMPERATURE')

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

# cd A_smoke_sensor

#for quantity in "${quantities[@]}"
#do
#  echo Slicefile quantity: "$quantity"
#  echo python3 preprocess_fds.py -q "$quantity"
#  python3 preprocess_fds.py -q "$quantity"
#done

python3 preprocess_fds.py -q 'CO2'
python3 preprocess_fds.py -q 'CO'
python3 preprocess_fds.py -q 'HCN'
python3 preprocess_fds.py -q 'TEMPERATURE'
python3 preprocess_fds.py -q 'EXTINCTION COEFFICIENT'
