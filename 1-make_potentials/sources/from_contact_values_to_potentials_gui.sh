#! /bin/bash -xe

model_name="$1"
symmetry="$2"
is_ipc="$3"

ls

pushd sources
  g++ -std=c++11 printPotential.cpp main.cpp -o compute.out
  
  target="../target_${model_name}_${symmetry}_contact"
  [ -d $target ] && rm -rf $target
  mkdir -p $target
  
  [ $ipc_model -eq 1 ] && is_ipc="-p"
  ./compute.out -c $is_ipc -m $symmetry -i inputfile -o $target
  mv inputfile ${target}/inputfile.dat
popd
