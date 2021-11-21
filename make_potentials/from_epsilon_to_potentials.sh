#!/bin/bash

# USAGE
# adapt the following variables to the model you want to simulate.
# do not delete any values, just substitute. keep the format

# name of the model you want to create (it is only for you)
model_name="45n"
# type of the model: can be janus, IPC (symmetric), or aIPC (asymmetric)
model="janus"  # janus    -> janus
               # ipc      -> symmetric IPC
               # asym_ipc -> asymmetric IPC
# delta (distance from the Hard Core at which the potential goes to zero)
delta=0.2
# patch eccentricity
ecc1=0.22
# epsilons
epsEE=0.245728
epsEP1=-3.11694
epsP1P1=21.2298
# potential normalization
emin=0.142495

# the following values are IGNORED for Janus and symmetric IPC
# patch_2 eccentricity
ecc2=0.22
epsEP2=-3.11694
epsP1P2=21.2298
epsP2P2=21.2298







# what follows is the implementation. you don't need to change it :)


pushd sources
  # generate inputfile
  echo $model_name > inputfile
  echo $delta >> inputfile
  echo $ecc1 >> inputfile
  echo $epsEE >> inputfile
  echo $epsEP1 >> inputfile
  echo $epsP1P1 >> inputfile
  echo $emin >> inputfile
  echo $ecc2 >> inputfile
  echo $epsEP2 >> inputfile
  echo $epsP1P2 >> inputfile
  echo $epsP2P2 >> inputfile
  g++ printPotential.cpp main.cpp -o compute.out

  mkdir -p ../target
  rm -rf ../target/lammpspot_${model_name}*

  ./compute.out -m $model -i inputfile -o ../target/lammpspot_${model_name}
popd 
