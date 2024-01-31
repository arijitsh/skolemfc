# Baseline Exact Counter

## Installing Dependencies

```
git submodule update --init

mkdir bin
cd deps

cd cryptominisat
mkdir build && cd build
cmake ..
make -j20
ln -s ./cryptominisat ../../bin/
cd ..

cd GPMC
./build.sh r
ln -s ./build/GPMC ../../bin/
```

## Running

```
./skolemcounter.sh myfile.qdimacs
```