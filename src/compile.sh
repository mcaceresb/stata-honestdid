#!/usr/bin/env bash

CWD=$(pwd -P)
ROOT=$(cd "$(dirname "$0")"; pwd -P)/..

[ ! -x "$(command -v git)"   ] && { echo "git executable not found";   exit 1; }
[ ! -x "$(command -v make)"  ] && { echo "make executable not found";  exit 1; }
[ ! -x "$(command -v cmake)" ] && { echo "cmake executable not found"; exit 1; }

cd ${ROOT}
git clone --recursive https://github.com/embotech/ecos
git clone --recursive https://github.com/osqp/osqp
cd ecos
make
cd ../osqp
mkdir build
cd build
cmake -G "Unix Makefiles" ..
cmake --build .
cd ../..
make all OSQP_H=./osqp/include OSQP_A=./osqp/build/out/libosqp.a ECOS_H="./ecos/include -I./ecos/external/SuiteSparse_config" ECOS_A="./ecos/libecos.a ./ecos/libecos_bb.a"
rm -rf ecos osqp
cd ${CWD}
exit 0
