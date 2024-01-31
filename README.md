# SkolemFC: An Approximate Skolem Function Counter

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- ![build](https://github.com/meelgroup/SkolemFC/workflows/build/badge.svg)
[![Docker Hub](https://img.shields.io/badge/docker-latest-blue.svg)](https://hub.docker.com/r/msoos/SkolemFC/) -->

SkolemFC is an approximate Skolem function counter.


## How to Build a Binary
To build on Linux, you will need the following:
```
sudo apt install build-essential cmake zlib1g-dev libboost-program-options-dev libboost-serialization-dev libgmp3-dev
```

Now clone this repository and run `./install.sh`, this should compile SkolemFC and all its dependencies.

```
git clone https://github.com/arijitsh/skolemfc/
cd skolemfc
./install.sh
```

Please follow `INSTALL.md` if the script reports some error, or you need more instructions for compiling in other OS, etc.



## How to Use the Binary
First, you must translate your problem to QDIMACS and just pass your file as input to SkolemFC, it will print the number of funcitons satisfying of the given QDIMACS formula.

### Running SkolemFC


```
$ ./skolemfc myfile.qdimacs

c [sklfc] SkolemFC Version: 36cf66a9ae
c [sklfc] executed with command line: ./skolemfc myfile.qdimacs
c [sklfc] using epsilon: 0.8 delta: 0.4 seed: 0
...
s fc 2 ** 4.00
...
c [sklfc] finished T: 0.25
c [sklfc] iterations: 729

```
SkolemFC reports that we have approximately `16 (=2 ** 4)` functions satisfying the QDIMACS specification.

### Guarantees
SkolemFC provides so-called "PAC", or Probably Approximately Correct, guarantees. In less fancy words, the system guarantees that the solution found is within a certain tolerance (called "epsilon") with a certain probability (called "delta"). The default tolerance and probability, i.e. epsilon and delta values, are set to 0.8 and 0.4, respectively. Both values are configurable.


### Issues, questions, bugs, etc.
Please click on "issues" at the top and [create a new issue](https://github.com/meelgroup/skolemfc/issues/new). All issues are responded to promptly.

## How to Cite

This work is by Arijit Shaw, Brendan Juba, and Kuldeep S. Meel, as [published in AAAI-24](https://arxiv.org/abs/2312.12026).

The benchmarks used in our evaluation can be found [here](https://zenodo.org/records/10404174).

## Exact Counter
An exact counter (termed "Baseline" in the paper) for Skolem Functions is available in the folder [`utils/baseline`](https://github.com/meelgroup/skolemfc/tree/main/utils/baseline). Please follow instructions in the README inside that folder for installing tools for that.
<!-- The old version, is available under the branch "paper". Please read the README of the old release to know how to compile the code. Old releases should easily compile. -->