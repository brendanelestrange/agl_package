#!/bin/bash
mkdir build
cd build 
rm -rf *
cmake ..
make 
mv agl_cpp.cpython-313-darwin.so ../src/ && \
cd ../src
