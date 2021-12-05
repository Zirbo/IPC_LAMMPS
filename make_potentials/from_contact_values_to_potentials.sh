#!/bin/bash

# USAGE
# adapt the following variables to the model you want to simulate.
# do not delete any values, just substitute. keep the format

# name of the model you want to create (it is only for you)
model_name="45n"
# type of the model: can be janus, IPC (symmetric), or aIPC (asymmetric)
model="ipc"
               # janus    -> janus
               # ipc      -> symmetric IPC
               # asym-ipc -> asymmetric IPC
# delta (distance from the Hard Core at which the potential goes to zero)
delta=0.2
# patch eccentricity
ecc1=0.22
# epsilons
vEE=0.1
vEP1=-1.0
vP1P1=4.0

# the following values are IGNORED for Janus and symmetric IPC
# patch_2 eccentricity
ecc2=0.22
vEP2=-1.0
vP1P2=4.0
vP2P2=4.0







# what follows is the implementation. you don't need to change it :)


pushd sources
  # generate inputfile
  echo $model_name > inputfile
  echo $delta >> inputfile
  echo $ecc1 >> inputfile
  echo $vEE >> inputfile
  echo $vEP1 >> inputfile
  echo $vP1P1 >> inputfile
  echo $ecc2 >> inputfile
  echo $vEP2 >> inputfile
  echo $vP1P2 >> inputfile
  echo $vP2P2 >> inputfile
  g++ printPotential.cpp main.cpp -o compute.out

  mkdir -p ../target
  rm -rf ../target/lammpspot_${model_name}*

  ./compute.out -c -m $model -i inputfile -o ../target/lammpspot_${model_name}
  mv inputfile ../target/inputfile_${model_name}_contact.dat
popd 
