#!/bin/sh

../configure --prefix=/usr/local/chroma/parscalar-gm-double --with-qdp=/usr/local/qdp++/parscalar-gm-double CXXFLAGS="" --enable-sse-wilson-dslash  --enable-opt-cfz-linop --enable-gmp  LIBS="-lgmp"
