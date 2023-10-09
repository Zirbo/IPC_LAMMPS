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

# parameters for the first patch


# delta -> TWICE the distance from the colloid surface (radius set to 0.5)
#               where the potential goes to zero
delta=0.2
# patch eccentricity
ecc1=0.22
# patch radius
rad1=0.38
# epsilons
epsEE=0.245728
epsEP1=-3.11694
epsP1P1=21.2298
# potential normalization
emin=0.142495

# parameters for the second patch -- **IGNORED** for Janus and symmetric
# patch_2 eccentricity
ecc2=0.22
# patch_2 radius
rad2=0.38
# other epsilons
epsEP2=-3.11694
epsP1P2=21.2298
epsP2P2=21.2298







# what follows is the implementation, you shouldn't change it


pushd sources
  # generate inputfile
  echo $model_name > inputfile
  echo $delta >> inputfile
  echo $ecc1 >> inputfile
  echo $rad1 >> inputfile
  echo $epsEE >> inputfile
  echo $epsEP1 >> inputfile
  echo $epsP1P1 >> inputfile
  echo $emin >> inputfile

  echo $ecc2 >> inputfile
  echo $rad2 >> inputfile
  echo $epsEP2 >> inputfile
  echo $epsP1P2 >> inputfile
  echo $epsP2P2 >> inputfile

  ./build.sh

  target="../target_${model_name}_${symmetry}_${mapping}_epsilons"
  [ -d $target ] && rm -rf $target
  mkdir -p $target

  [ $ipc_model -eq 1 ] && is_ipc="-p"
  ./bld/lammps_pot_generator \
    -e $is_ipc -s $symmetry -m $mapping \
    -i inputfile -o $target
  rm inputfile
popd 
