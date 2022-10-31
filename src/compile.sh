#!/usr/bin/env bash

CWD=$(pwd -P)
ROOT=$(cd "$(dirname "$0")"; pwd -P)/..

[ ! -x "$(command -v git)"   ] && { echo "git executable not found";   exit 1; }
[ ! -x "$(command -v make)"  ] && { echo "make executable not found";  exit 1; }
[ ! -x "$(command -v cmake)" ] && { echo "cmake executable not found"; exit 1; }
[ ! -x "$(command -v sed)"   ] && { echo "sed executable not found";   exit 1; }

if [[ $(uname -s) == *'CYGWIN'* ]]; then
    [ ! -x "$(command -v x86_64-w64-mingw32-gcc)"   ] && { echo "x86_64-w64-mingw32-gcc not found";   exit 1; }
fi

cd ${ROOT}
git clone --recursive https://github.com/embotech/ecos
git clone --recursive https://github.com/osqp/osqp
cd ecos
if [[ $(uname -s) == *'CYGWIN'* ]]; then
    make CC=x86_64-w64-mingw32-gcc
else
    make
fi
# [ $(uname -s) == 'Darwin' ] && sed -i -e 's/^CFLAGS\(.*\)/CFLAGS\1 -arch arm64/g' Makefile
cd ../osqp
if [[ $(uname -m) == 'arm64' && $(uname -s) == 'Darwin' ]]; then
    sed -i -e '/^elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin").*/a\'$'\n''set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE INTERNAL "" FORCE)\
        ' ../CMakeLists.txt
fi
mkdir build
cd build
if [[ $(uname -s) == *'CYGWIN'* ]]; then
    CC=x86_64-w64-mingw32-gcc cmake -G "Unix Makefiles" -DIS_WINDOWS=ON ..
    CC=x86_64-w64-mingw32-gcc cmake --build .
else
    cmake -G "Unix Makefiles" ..
    cmake --build .
fi
cd ../..
if [[ $(uname -m) == 'arm64' && $(uname -s) == 'Darwin' ]]; then
    make all OSQP_H=./osqp/include OSQP_A=./osqp/build/out/libosqp.a ECOS_H="./ecos/include -I./ecos/external/SuiteSparse_config" ECOS_A="./ecos/libecos.a ./ecos/libecos_bb.a" OSQP_OUT=src/build/honestosqp_macosxarm64.plugin ECOS_OUT=src/build/honestecos_macosxarm64.plugin
else
    make all OSQP_H=./osqp/include OSQP_A=./osqp/build/out/libosqp.a ECOS_H="./ecos/include -I./ecos/external/SuiteSparse_config" ECOS_A="./ecos/libecos.a ./ecos/libecos_bb.a"
fi
rm -rf ecos osqp
cd ${CWD}
exit 0
