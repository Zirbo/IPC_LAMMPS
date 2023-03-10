#!/bin/bash

mkdir -p bld
pushd bld
  cmake ../sources/
  make -j4
popd
cp bld/lammpsIPCpostprocess lammpsIPCpostprocess.out
