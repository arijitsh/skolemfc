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
ln -s $(pwd)/cryptominisat5 ../../../bin/
cd ../..

cd GPMC
./build.sh r
ln -s $(pwd)/bin/gpmc ../../bin/
cd ..
```

## Running

```
./skolemcounter.sh myfile.qdimacs
```