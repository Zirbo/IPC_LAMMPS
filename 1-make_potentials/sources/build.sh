#! /bin/bash

eigen=$( ls ../../dependencies/eigen/ )

if [ -z "$eigen" ]; then
  pushd ../..
    git submodule init
    git submodule update
  popd
fi

mkdir -p bld
pushd bld
  cmake ../
  make -j4
popd
