#!/usr/bin/bash

git submodule update --init --recursive
cd deps

cd louvain-community
mkdir build && cd build
cmake ..
make -j20
cd ..

cd cryptominisat
mkdir build && cd build
cmake ..
make -j20
cd ..

cd arjun
mkdir build && cd build
cmake ..
make -j20
cd ..

cd approxmc
mkdir build && cd build
cmake ..
make -j20
cd ..

cd unigen
mkdir build && cd build
cmake ..
make -j20
cd ../..

cd ..
mkdir build && cd build
cmake ..
make -j20

