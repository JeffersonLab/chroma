#!/bin/sh

../configure --prefix=$HOME/arch/chroma/scalar --with-qdp=/usr/local/qdp++/qdp1-20-5/scalar CXXFLAGS="" --enable-sse-wilson-dslash --enable-opt-cfz-linop --with-gmp=/usr  --enable-testcase-runner=trivial

