#!/bin/bash

set -e

command rm -rf bin build lib 
mkdir build
cd build
#FC=gfortran cmake .. -DCMAKE_BUILD_TYPE=DEBUG 
#-DLAPACK_LIBRARIES= -DBLAS_LIBRARIES=
FC=ifort cmake .. -DCMAKE_BUILD_TYPE=DEBUG
#FC=ifort cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make
