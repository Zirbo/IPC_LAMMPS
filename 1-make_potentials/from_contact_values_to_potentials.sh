#!/bin/bash

# USAGE
# adapt the following variables to the model you want to simulate.
# do not delete any values, just substitute. keep the format

# name of the model you want to create
# it will only appear on the output directory, so choose freely
model_name="your-model-name"

# type of the model you want to create.
# janus models have a single patch, symmetric models a second symmetric one
# asymmetric models have two different patches, so you will need to specify
# the parameters of the second one below
# recognized values are:
# j-ospc -> janus colloid
# s-ospc -> two-symmetric-patches colloid
# a-ospc -> two-aymmetric-patches colloid
# j-ipc -> janus IPC
# s-ipc -> symmetric IPC
# a-ipc -> asymmetric IPC
ospc_type="s-ipc"

# parameters for the first patch

# delta (distance from the Hard Core at which the potential goes to zero)
delta=0.2
# patch eccentricity
ecc1=0.22
# patch radius
rad1=0.38
# contact values
vEE=0.1
vEP1=-1.0
vP1P1=4.0

#parameters for the second patch -- IGNORED for Janus and symmetric
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

  target="../target_${model_name}_${ospc_type}_contact"
  [ -d $target ] && rm -rf $target
  mkdir -p $target

  ./compute.out -c -m $ospc_type -i inputfile -o $target
  mv inputfile ${target}/inputfile.dat
popd 
