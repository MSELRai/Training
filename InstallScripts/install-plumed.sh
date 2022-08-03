#!/bin/bash

#./configure --prefix=$(pwd)/install --disable-mpi

version=2.8.0

wget https://github.com/plumed/plumed2/releases/download/v${version}/plumed-${version}.tgz
tar -xzvf plumed-${version}.tgz
rm plumed-${version}.tgz
mv plumed-${version} ${version}

cd $version

OPT="-O2 -fPIC"
./configure --prefix=$(pwd)/install \
            CC=mpicc \
            CXX=mpic++ \
            CFLAGS="${OPT}" \
            CXXFLAGS="${OPT}" \
            MPIEXEC="mpirun" \
            PYTHON_BIN=$(which python3)

make -j2 && make install
