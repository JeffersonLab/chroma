#!/bin/sh

../configure --prefix=/usr/local/chroma/scalar --with-qdp=/usr/local/qdp++/scalar CXXFLAGS="" --enable-sse-wilson-dslash --enable-gmp LIBS="-lgmp"

