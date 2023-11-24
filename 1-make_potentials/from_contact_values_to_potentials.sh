#!/bin/bash

# USAGE
# adapt the following variables to the model you want to simulate.
# do not delete any values, just substitute. keep the format

# name of the model you want to create
# it will only appear on the output directory, so choose freely
model_name="your-model-name"

# type of mapping, two are available:
# g -> geometrical, from the original paper
# e -> exponential, from the paper with ....
mapping=g

# symmetry of the colloid that you want to create.
#              accepted values:
# janus -> janus colloids with a single patch
# symm -> colloids with two polar and identical patches
# asymm -> colloids with two polar but different (in geometry and/or interaction) patches
#        for asymm you must also specify the parameters of the second patch
symmetry="asymm"

# if you want to simulate an ipc model, set this value to 1
# any other value disables it
ipc_model=0

# lower threshold value for both potential and force
cutoff=0.0001

# parameters for the first patch


# delta -> TWICE the distance from the colloid surface (radius set to 0.5)
#               where the potential goes to zero
delta=0.2
# patch eccentricity
ecc1=0.22
# patch radius
rad1=0.38

# contact values -> the naming is referring to two patch systems.
vEE=0.1     # for JANUS -> backside-backside
vEP1=-1.0   # for JANUS -> backside-patch
vP1P1=4.0   # for JANUS -> patch-patch




# parameters for the second patch -- **IGNORED** for Janus and symmetric
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
  echo $cutoff >> inputfile
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

  ./build.sh

  target="../target_${model_name}_${symmetry}_${mapping}_contact"
  [ -d $target ] && rm -rf $target
  mkdir -p $target

  [ $ipc_model -eq 1 ] && is_ipc="-p"
  ./bld/lammps_pot_generator \
    -c $is_ipc -s $symmetry -m $mapping \
    -i inputfile -o $target
  rm inputfile
popd 
