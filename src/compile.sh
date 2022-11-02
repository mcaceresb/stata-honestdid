#!/usr/bin/env bash

CWD=$(pwd -P)
ROOT=$(cd "$(dirname "$0")"; pwd -P)/..

[ ! -x "$(command -v git)"   ] && { echo "git executable not found; please install git";     exit 1; }
[ ! -x "$(command -v make)"  ] && { echo "make executable not found; please install make";   exit 1; }
[ ! -x "$(command -v cmake)" ] && { echo "cmake executable not found; please install cmake"; exit 1; }
[ ! -x "$(command -v sed)"   ] && { echo "sed executable not found; please install sed";     exit 1; }

if [[ $(uname -s) == *'CYGWIN'* ]]; then
    DETECT_OS=CYGWIN
elif [[ $(uname -s) == 'Darwin' && $(uname -m) == 'arm64' ]]; then
    DETECT_OS=APPLE_ARM64
elif [[ $(uname -s) == 'Darwin' && $(uname -m) == 'x86_64' ]]; then
    DETECT_OS=APPLE_X86_64
elif [[ $(uname -s) == 'Linux' ]]; then
    DETECT_OS=LINUX
else
    DETECT_OS=UNKNOWN
fi
: RUN_OS=${RUN_OS:=${DETECT_OS}}
: RUN_CLONE=${RUN_CLONE:=1}
: RUN_BUILD=${RUN_BUILD:=1}
: RUN_COPY=${RUN_COPY:=1}

if [[ "${RUN_OS}" != "${DETECT_OS}" ]]; then
    echo "OS detected as '${DETECT_OS}' but running as ${RUN_OS}"
elif [[ "${RUN_OS}" == "UNKNOWN" ]]; then
    echo "Unknown OS $(uname -s)"
    exit 1
else
    echo "Running as ${RUN_OS} ($(uname -s) $(uname -m))"
fi

[[ "${RUN_OS}" == "CYGWIN" && ! -x "$(command -v x86_64-w64-mingw32-gcc)" ]] && { echo "x86_64-w64-mingw32-gcc not found"; exit 1; }

cd ${ROOT}
(( ${RUN_CLONE} )) && git clone --recursive https://github.com/embotech/ecos
(( ${RUN_CLONE} )) && git clone --recursive https://github.com/osqp/osqp

echo "Compiling ECOS plugin"
cd ecos
make purge
cp ecos.mk ecos.mk.bak
[[ "${RUN_OS}" == "APPLE_ARM64"  ]] && sed -i -e 's/^CFLAGS\(.*\)/CFLAGS\1 -arch arm64/'     ecos.mk
[[ "${RUN_OS}" == "APPLE_X86_64" ]] && sed -i -e 's/^CFLAGS\(.*\)/CFLAGS\1 -arch x86_64/'    ecos.mk
[[ "${RUN_OS}" == *"APPLE"*      ]] && sed -i -e 's/^#CC = .*/CC = clang/'                   ecos.mk
[[ "${RUN_OS}" == "CYGWIN"       ]] && sed -i -e 's/^#CC = .*/CC = x86_64-w64-mingw32-gcc/'  ecos.mk
[[ "${RUN_OS}" == "CYGWIN"       ]] && sed -i -e 's/^ISWINDOWS := 0/ISWINDOWS := 1/g'        ecos.mk
(( ${RUN_BUILD} )) && make
cp -f ecos.mk.bak ecos.mk
rm -f ecos.mk.bak

echo "Compiling OSQP plugin"
cd ../osqp
cp CMakeLists.txt CMakeLists.txt.bak
CMAKE_FLAGS=
if [[ "${RUN_OS}" == "APPLE_X86_64" ]]; then
    sed -i -e '/^elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin").*/a\'$'\n''set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE INTERNAL "" FORCE)\
        ' CMakeLists.txt
elif [[ "${RUN_OS}" == "APPLE_ARM64" ]]; then
    sed -i -e '/^elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin").*/a\'$'\n''set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE INTERNAL "" FORCE)\
        ' CMakeLists.txt
elif [[ "${RUN_OS}" == "CYGWIN" ]]; then
    CMAKE_FLAGS="-DIS_WINDOWS=ON"
fi
[[ "${RUN_OS}" == *"APPLE"* ]] && echo "set(CMAKE_C_COMPILER "clang")"                  >>  CMakeLists.txt
[[ "${RUN_OS}" == "CYGWIN"  ]] && echo "set(CMAKE_C_COMPILER "x86_64-w64-mingw32-gcc")" >>  CMakeLists.txt

[[ "${CMAKE_FLAGS}" != "" ]] && echo "running cmake with added flags ${CMAKE_FLAGS}"
[ -d "build" ] && rm -rf build
mkdir build && cd build
cmake ${CMAKE_FLAGS} -G "Unix Makefiles" ..
CMAKE_FLAGS=
(( ${RUN_BUILD} )) && cmake --build .
cd ..
cp -f CMakeLists.txt.bak CMakeLists.txt
rm -f CMakeLists.txt.bak
cd ..

if [[ "${RUN_OS}" == "APPLE_X86_64" ]]; then
    HONEST_OUT="OSQP_OUT=src/build/honestosqp_macosx86_64.plugin ECOS_OUT=src/build/honestecos_macosx86_64.plugin GCC=clang"
    HONEST_FLAGS="-bundle -DSYSTEM=APPLEMAC -arch x86_64"
elif [[ "${RUN_OS}" == "APPLE_ARM64"  ]]; then
    HONEST_OUT="OSQP_OUT=src/build/honestosqp_macosxarm64.plugin ECOS_OUT=src/build/honestecos_macosxarm64.plugin GCC=clang"
    HONEST_FLAGS="-bundle -DSYSTEM=APPLEMAC -arch arm64"
elif [[ "${RUN_OS}" == "CYGWIN"  ]]; then
    HONEST_OUT=
    HONEST_FLAGS="-shared -fPIC"
elif [[ "${RUN_OS}" == "LINUX"  ]]; then
    HONEST_OUT=
    HONEST_FLAGS="-shared -fPIC -DSYSTEM=OPUNIX"
fi

echo "Compiling honestdid plugin"
[[ "${HONEST_FLAGS}" != "" ]] && echo "    with flags ${HONEST_FLAGS}"
[[ "${HONEST_OUT}"   != "" ]] && echo "    with opts  ${HONEST_OUT}"
if (( ${RUN_BUILD} )); then
    make all OSQP_H=./osqp/include OSQP_A=./osqp/build/out/libosqp.a ECOS_H="./ecos/include -I./ecos/external/SuiteSparse_config" ECOS_A="./ecos/libecos.a ./ecos/libecos_bb.a" ${HONEST_OUT} OSFLAGS="${HONEST_FLAGS}"
    HONEST_RC=$?
    [[ ${HONEST_RC} -eq 0 ]] && echo "Plugin compiled successfully" || echo "Plugin compiletion failed; please report issue"
else
    HONEST_RC=0
fi

if (( ${RUN_COPY} )) && [[ "${RUN_OS}" == *"APPLE"* ]]; then
    echo "Overriding OSX plugin with locally compiled plugin"
    eval `echo ${HONEST_OUT}`
    cp -f ${OSQP_OUT} src/build/honestosqp_macosx.plugin
    cp -f ${ECOS_OUT} src/build/honestecos_macosx.plugin
fi

(( ${RUN_CLONE} )) && rm -rf ecos osqp
cd ${CWD}
exit ${HONEST_RC}
