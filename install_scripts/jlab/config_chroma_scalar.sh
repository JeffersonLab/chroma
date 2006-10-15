#!/bin/sh

../configure --prefix=$HOME/arch/chroma/scalar --with-qdp=/dist/scidac/qdp++/qdp1-21-3/scalar CXXFLAGS="" --enable-sse-wilson-dslash --enable-opt-cfz-linop --with-gmp=/usr  --enable-testcase-runner=trivial

