#!/bin/bash

set -e

command rm -rf bin build lib 
mkdir build
cd build
#FC=gfortran cmake .. -DCMAKE_BUILD_TYPE=DEBUG 
#-DLAPACK_LIBRARIES= -DBLAS_LIBRARIES=
FC=ifx cmake .. -DCMAKE_BUILD_TYPE=DEBUG
#FC=ifx cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make
