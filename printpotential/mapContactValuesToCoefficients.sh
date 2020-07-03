#!/bin/bash

# name of the model you want to create
model_name=45n
# delta (distance from the Hard Core at which the potential goes to zero)
delta=0.2
# patch eccentricity
ecc=0.22
# potential at contact in the equator-equator alignment
vEE=0.1
# potential at contact in the patch-equator alignment
vEP=-1.0
# potential at contact in the patch-patch alignment
vPP=4.0




pushd sources
  python3 mapContactValuesToCoefficients.py $model_name $delta $ecc $vEE $vEP $vPP
  if [ ! -f compute.out ]; then
    g++ printPotential.cpp -o compute.out
  fi
  ./compute.out inputfile_${model_name}.txt ../lammpspot_${model_name}
  mv inputfile_${model_name}.txt ..
  mv MC_inputfile_${model_name}.txt ..
popd 

# ADVANCED
# if you want to change other parameters than the one specified here,
# modify your inputfile_${model_name}.txt (a commented example is in
# sources/inputfile_example.txt) and rerun
# ../compute.out inputfile_${model_name}.txt lammpspot_${model_name}
