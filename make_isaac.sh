#!/bin/bash
module load eigen/3.4.0-gcc
module load cmake/3.30.5-gcc
mkdir build
cd build 
rm -rf *
cmake ..
make 
mv *.so ../src/
cd ../
