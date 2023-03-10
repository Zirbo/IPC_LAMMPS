#! /bin/bash

eigen=$( ls ../../dependencies/eigen/ )

if [ -z "$eigen" ]; then
  pushd ../..
    git submodule init
    git submodule update
  popd
fi

g++ -std=c++14 -I ../../dependencies/eigen/ printPotential.cpp main.cpp -o compute.out
