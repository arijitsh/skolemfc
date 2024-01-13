# SkolemFC: An Approximate Skolem Function Counter

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- ![build](https://github.com/meelgroup/SkolemFC/workflows/build/badge.svg)
[![Docker Hub](https://img.shields.io/badge/docker-latest-blue.svg)](https://hub.docker.com/r/msoos/SkolemFC/) -->

SkolemFC is an approximate Skolem function counter.


## How to Build a Binary
To build on Linux, you will need the following:
```
sudo apt-get install build-essential cmake
sudo apt-get install zlib1g-dev libboost-program-options-dev libboost-serialization-dev
apt-get install libgmp3-dev
```

Then, build CryptoMiniSat, Arjun, and SkolemFC:
```
git clone https://github.com/msoos/cryptominisat
cd cryptominisat
mkdir build && cd build
cmake ..
make
sudo make install
sudo ldconfig

cd ../..
git clone https://github.com/meelgroup/arjun
cd arjun
mkdir build && cd build
cmake ..
make
sudo make install
sudo ldconfig

cd ../..
git clone https://github.com/meelgroup/approxmc
cd approxmc
mkdir build && cd build
cmake ..
make
sudo make install
sudo ldconfig


cd ../..
git clone https://github.com/meelgroup/skolemfc
cd skolemfc
mkdir build && cd build
cmake ..
make
sudo make install
sudo ldconfig
```


## How to Use the Binary
First, you must translate your problem to QDIMACS and just pass your file as input to SkolemFC, it will print the number of funcitons satisfying of the given QDIMACS formula.

### Running SkolemFC


```
$ SkolemFC --seed 5 myfile.cnf
c SkolemFC version 0.0.1
[appmc] FINISHED SkolemFC T: 0.04 s
c [appmc] Number of solutions is: 2 ** 4
s fc 2 ** 4
```
SkolemFC reports that we have approximately `16 (=2 ** 4)` functions satisfying the QDIMACS specification.

### Guarantees
SkolemFC provides so-called "PAC", or Probably Approximately Correct, guarantees. In less fancy words, the system guarantees that the solution found is within a certain tolerance (called "epsilon") with a certain probability (called "delta"). The default tolerance and probability, i.e. epsilon and delta values, are set to 0.8 and 0.2, respectively. Both values are configurable.


### Issues, questions, bugs, etc.
Please click on "issues" at the top and [create a new issue](https://github.com/meelgroup/skolemfc/issues/new). All issues are responded to promptly.

## How to Cite

This work is by Arijit Shaw, Brendan Juba, and Kuldeep S. Meel, as [published in AAAI-24](https://arxiv.org/abs/2312.12026).

The benchmarks used in our evaluation can be found [here](https://zenodo.org/records/10404174).

## Old Versions
The old version, is available under the branch "paper". Please read the README of the old release to know how to compile the code. Old releases should easily compile.