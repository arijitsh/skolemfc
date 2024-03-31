#!/usr/bin/bash
set -x
set -e

git submodule update --init
(
cd deps

(cd louvain-community
mkdir -p build && cd build
cmake -DSTATICCOMPILE=ON ..
make -j20
cd ..)

(cd cryptominisat
mkdir -p build && cd build
cmake -DSTATICCOMPILE=ON ..
make -j20
cd ..)

(cd arjun
mkdir -p build && cd build
cmake -DSTATICCOMPILE=ON ..
make -j20
cd ..)

(cd approxmc
mkdir -p build && cd build
cmake -DSTATICCOMPILE=ON ..
make -j20
cd ..)

(cd unigen
mkdir -p build && cd build
cmake -DSTATICCOMPILE=ON ..
make -j20
cd ..)

cd ..
mkdir -p build && cd build
cmake -DSTATICCOMPILE=ON ..
make -j20

)
