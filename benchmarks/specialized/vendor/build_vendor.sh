#!/bin/sh
# Build GraSS and its bsc backend inside the locked image.
# Invoked by fetch_vendor.ps1; assumes the project is mounted at /work.
set -eu

V=/work/benchmarks/specialized/vendor
mkdir -p "$V/bin" "$V/libbsc/build" "$V/GraSS/build"

cd "$V/libbsc/build"
g++ -O2 -c \
    ../libbsc/adler32/adler32.cpp \
    ../libbsc/bwt/libsais/libsais.c \
    ../libbsc/bwt/bwt.cpp \
    ../libbsc/coder/coder.cpp \
    ../libbsc/coder/qlfc/qlfc.cpp \
    ../libbsc/coder/qlfc/qlfc_model.cpp \
    ../libbsc/filters/detectors.cpp \
    ../libbsc/filters/preprocessing.cpp \
    ../libbsc/libbsc/libbsc.cpp \
    ../libbsc/lzp/lzp.cpp \
    ../libbsc/platform/platform.cpp \
    ../libbsc/st/st.cpp \
    ../bsc.cpp
g++ -o bsc \
    bsc.o adler32.o bwt.o coder.o detectors.o libbsc.o libsais.o \
    lzp.o platform.o preprocessing.o qlfc.o qlfc_model.o st.o
cp bsc "$V/bin/bsc"

cd "$V/GraSS"
javac -d build src/*.java
