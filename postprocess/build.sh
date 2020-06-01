#!/bin/bash

mkdir -p bld
pushd bld
  cmake ..
  make -j4
popd
