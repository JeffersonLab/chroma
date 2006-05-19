#!/bin/tcsh

../configure --prefix=/usr/local/share/chroma/scalar --with-qdp=/usr/local/share/qdp++/scalar CFLAGS="-O2 -msse -msse2 -march=pentium4" CXXFLAGS="" LIBS="-lgmp" --enable-sse-wilson-dslash  --with-gmp=/usr --enable-opt-cfz-linop  --enable-testcase-runner=trivial

