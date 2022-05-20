#!/bin/bash

mkdir -p bld
pushd bld
  cmake ../sources/
  make -j4
popd
