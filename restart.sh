#!/bin/bash

set -e

command rm -rf bin build lib 
mkdir build
cd build
FC=ifort cmake .. -DCMAKE_BUILD_TYPE=DEBUG
#FC=ifort cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make
