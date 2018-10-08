#!/bin/sh

# Execution of preprocessing script to supply JPSfire

#===================CHANGE PARAMETERS AS NECEASSARY=============================

quantities_FED=(
'CARBON DIOXIDE VOLUME FRACTION'
'CARBON MONOXIDE VOLUME FRACTION'
'HYDROGEN CYANIDE VOLUME FRACTION'
'TEMPERATURE'
)

locations=(
'Z',2.25 #height for extinction coefficient (must be compliant with FDS simulation)
'Z',1.75 #height for FED calculation (must be compliant with FDS simulation)
)
#==================AUTOMATIC PART - NO CHANGES NECEASSARY=======================

#FED grids
IFS=","; set ${locations[1]};
# $1 = slice_dim as str and $2 = slice_coord as float

for quantity in "${quantities_FED[@]}"
do
  echo Slicefile quantity: "$quantity"
  echo Slicefile dimension: "$1", coordinate: $2 m
  python3 preprocess_fds.py -q "$quantity" -d "$1" -c $2 -p
done

#extinction grids
python3 preprocess_fds.py -p
