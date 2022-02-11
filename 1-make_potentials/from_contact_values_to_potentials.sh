#!/bin/bash

# USAGE
# adapt the following variables to the model you want to simulate.
# do not delete any values, just substitute. keep the format

# name of the model you want to create
# it will only appear on the output directory, so choose freely
model_name="your-model-name"

# symmetry of the colloid that you want to create.
#              accepted values:
# janus -> janus colloids with a single patch
# symm -> symmetric symmetric colloids with two symmetric patches
# asymm -> asymmetric colloids with two different patches
#        for asymm you must also specify the parameters of the second patch
symmetry="asymm"

# if you want to simulate an ipc model, set this value to 1
# any other value disables
ipc_model=0

# parameters for the first patch

# delta (twice the distance from the Hard Core where the potential goes to zero)
delta=0.2
# patch eccentricity
ecc1=0.22
# patch radius
rad1=0.38
# contact values
vEE=0.1
vEP1=-1.0
vP1P1=4.0

#parameters for the second patch -- **IGNORED** for Janus and symmetric
# patch_2 eccentricity
ecc2=0.22
# patch_2 radius
rad2=0.38
# other contact values
vEP2=-1.0
vP1P2=4.0
vP2P2=4.0







# what follows is the implementation, you shouldn't change it


pushd sources
  # generate inputfile
  echo $model_name > inputfile
  echo $delta >> inputfile
  echo $ecc1 >> inputfile
  echo $rad1 >> inputfile
  echo $vEE >> inputfile
  echo $vEP1 >> inputfile
  echo $vP1P1 >> inputfile

  echo $ecc2 >> inputfile
  echo $rad2 >> inputfile
  echo $vEP2 >> inputfile
  echo $vP1P2 >> inputfile
  echo $vP2P2 >> inputfile
  g++ -std=c++11 printPotential.cpp main.cpp -o compute.out

  target="../target_${model_name}_${symmetry}_contact"
  [ -d $target ] && rm -rf $target
  mkdir -p $target

  [ $ipc_model -eq 1 ] && is_ipc="-p"
  ./compute.out -c ${is_ipc} -m $symmetry -i inputfile -o $target
  mv inputfile ${target}/inputfile.dat
popd 
